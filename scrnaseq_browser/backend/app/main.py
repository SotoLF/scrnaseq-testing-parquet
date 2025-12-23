from fastapi import FastAPI, HTTPException, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import numpy as np
import hashlib
import base64
import pandas as pd

from .parquet_backend import get_backend
from .config import get_datasets
from .cache import get_json, set_json, get_compressed, set_compressed

app = FastAPI(title="scRNA Parquet API", version="2.0")

@app.on_event("startup")
async def startup_event():
    print("\n" + "="*40)
    print("       Backend Startup")
    print("="*40)
    datasets = get_datasets()
    print(f"Loaded {len(datasets)} datasets from configuration:")
    for d in datasets:
        print(f" - {d.id}: {d.name} ({d.uri})")
    print("="*40 + "\n")

# Enable CORS for Shiny app
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class ExprRequest(BaseModel):
    dataset_id: str
    measurement: str = "RNA"
    layer: str = "X"
    gene: str
    idx0: List[int]  # obs soma_joinid values (NOT row positions)

class BatchRequest(BaseModel):
    dataset_id: str
    measurement: str = "RNA"
    layer: str = "X"
    genes: List[str]
    idx0: List[int]  # obs soma_joinid values (NOT row positions)

def b64_f32(arr: np.ndarray) -> str:
    """Convert float32 array to base64 string"""
    return "base64:" + base64.b64encode(arr.astype(np.float32).tobytes()).decode('ascii')

def b64_i32(arr: np.ndarray) -> str:
    """Convert int32 array to base64 string"""
    return "base64:" + base64.b64encode(arr.astype(np.int32).tobytes()).decode('ascii')

def _hash_idx(idx0: List[int]) -> str:
    """Create hash of cell indices for cache keys"""
    b = np.asarray(idx0, dtype=np.int32).tobytes()
    return hashlib.md5(b).hexdigest()

@app.get("/")
def root():
    """API root endpoint"""
    try:
        # Just show info about the first dataset or general info
        datasets_list = get_datasets()
        return {
            "message": "scRNA-seq Parquet Backend API",
            "version": "2.0",
            "backend": "Parquet",
            "datasets": [d.id for d in datasets_list],
            "endpoints": [
                "/datasets",
                "/genes/{dataset_id}",
                "/obs_types/{dataset_id}",
                "/obsm_xy/{dataset_id}",
                "/obs_col/{dataset_id}/{column}",
                "/expr_subset",
                "/expr_batch"
            ]
        }
    except Exception as e:
        return {
            "message": "scRNA-seq Parquet Backend API",
            "version": "2.0",
            "error": str(e)
        }

@app.get("/datasets")
def datasets():
    """List all available datasets"""
    out = []
    for cfg in get_datasets():
        try:
            # Initialize backend to get counts (lazy load)
            backend = get_backend(cfg.id)
            out.append({
                "id": cfg.id,
                "name": cfg.name,
                "n_cells": backend.obs_count(),
                "n_genes": backend.var_count()
            })
        except Exception as e:
            print(f"Error loading dataset {cfg.id}: {e}")
            # Include it but with error info or skip?
            # Better to skip or show error state
            out.append({
                "id": cfg.id,
                "name": f"{cfg.name} (Error)",
                "n_cells": 0,
                "n_genes": 0
            })
    return out

@app.get("/info/{dataset_id}")
def info(dataset_id: str):
    """Get dataset metadata, layers, embeddings, etc."""
    backend = get_backend(dataset_id)
    # Match what the Shiny client expects.
    return {
        "id": dataset_id,
        "n_cells": backend.obs_count(),
        "n_genes": backend.var_count(),
        "measurement": "RNA",
        "layers": backend.available_layers,
        "embeddings": backend.available_embeddings
    }



@app.get("/genes/{dataset_id}")
def genes(dataset_id: str, measurement: str = "RNA"):
    """Get list of all genes"""
    key = f"genes:{dataset_id}:{measurement}"
    cached = get_json(key)
    if cached is not None:
        return cached
    backend = get_backend(dataset_id)
    out = backend.genes()
    set_json(key, out, ex=24*3600)
    return out

@app.get("/obs_types/{dataset_id}")
def obs_types(dataset_id: str):
    """Get obs column types"""
    key = f"obs_types:{dataset_id}"
    cached = get_json(key)
    if cached is not None:
        return cached
    backend = get_backend(dataset_id)
    out = backend.get_obs_types()
    set_json(key, out, ex=24*3600)
    return out

@app.get("/obsm_xy/{dataset_id}")
def obsm_xy(dataset_id: str, measurement: str = "RNA", key: Optional[str] = None):
    """
    Returns embedding coordinates with binary compressed caching.
    """
    if key is None:
        key = "X_umap"
    cache_key = f"obsm_xy:{dataset_id}:{measurement}:{key}"
    
    # Try binary cache first
    cached_bytes = get_compressed(cache_key)
    if cached_bytes is not None:
        arr = np.frombuffer(cached_bytes, dtype=np.float32).reshape(-1, 2)
        return {"x": b64_f32(arr[:, 0]), "y": b64_f32(arr[:, 1]), "n": int(arr.shape[0])}
    
    # Get data from Parquet backend
    backend = get_backend(dataset_id)
    x, y = backend.get_obsm_xy(key)
    arr = np.column_stack([x, y]).astype(np.float32)
    set_compressed(cache_key, arr.tobytes(), ex=24*3600)
    
    return {"x": b64_f32(arr[:, 0]), "y": b64_f32(arr[:, 1]), "n": int(arr.shape[0])}

@app.get("/obs_col/{dataset_id}/{col}")
def obs_col(dataset_id: str, col: str):
    """Get obs column data in a compact typed format expected by the Shiny client."""
    meta_key = f"obs_col:{dataset_id}:{col}:meta"
    data_key = f"obs_col:{dataset_id}:{col}:data"

    cached_meta = get_json(meta_key)
    cached_bytes = get_compressed(data_key)
    if cached_meta is not None and cached_bytes is not None:
        if cached_meta["type"] == "num":
            arr = np.frombuffer(cached_bytes, dtype=np.float32)
            return {"type": "num", "data": b64_f32(arr)}
        else:
            codes = np.frombuffer(cached_bytes, dtype=np.int32)
            return {
                "type": "cat",
                "codes": b64_i32(codes),
                "levels": cached_meta["levels"],
            }

    backend = get_backend(dataset_id)
    try:
        values = backend.get_obs_column(col)
    except Exception as e:
        raise HTTPException(status_code=404, detail=f"obs column not found: {col}") from e

    # Pandas categorical (fast path)
    if isinstance(values, pd.Categorical) or pd.api.types.is_categorical_dtype(getattr(values, "dtype", None)):
        codes = np.asarray(getattr(values, "codes"), dtype=np.int32)
        levels = [str(x) for x in list(getattr(values, "categories"))]
        set_compressed(data_key, codes.tobytes(), ex=6 * 3600)
        set_json(meta_key, {"type": "cat", "levels": levels}, ex=6 * 3600)
        return {"type": "cat", "codes": b64_i32(codes), "levels": levels}

    # Numeric
    dtype = getattr(values, "dtype", np.asarray(values).dtype)
    if pd.api.types.is_numeric_dtype(dtype) or pd.api.types.is_bool_dtype(dtype):
        arr = np.asarray(values, dtype=np.float32)
        set_compressed(data_key, arr.tobytes(), ex=6 * 3600)
        set_json(meta_key, {"type": "num"}, ex=6 * 3600)
        return {"type": "num", "data": b64_f32(arr)}

    # Categorical/string/object
    as_str = np.asarray(values).astype(str)
    levels, codes = np.unique(as_str, return_inverse=True)
    codes = codes.astype(np.int32)
    levels_list = levels.tolist()

    set_compressed(data_key, codes.tobytes(), ex=6 * 3600)
    set_json(meta_key, {"type": "cat", "levels": levels_list}, ex=6 * 3600)

    return {
        "type": "cat",
        "codes": b64_i32(codes),
        "levels": levels_list,
    }

@app.post("/expr_subset")
def expr_subset(req: ExprRequest):
    """
    FAST PARQUET: Get expression values for a gene.
    Now <100ms instead of 5-7 seconds!
    """
    h = _hash_idx(req.idx0)
    cache_key = f"expr:{req.dataset_id}:{req.measurement}:{req.layer}:{req.gene}:{h}"
    
    # Try binary cache first
    cached_bytes = get_compressed(cache_key)
    if cached_bytes is not None:
        arr = np.frombuffer(cached_bytes, dtype=np.float32)
        return {"data": b64_f32(arr)}
    
    # Get data from Parquet backend (instant!)
    backend = get_backend(req.dataset_id)
    idx = np.asarray(req.idx0, dtype=np.int64) if req.idx0 else None
    
    try:
        arr = backend.get_gene_expression(req.gene, cell_indices=idx, layer=req.layer)
    except ValueError:
        # Gene not found or similar
        return {"data": None}
    except Exception as e:
        print(f"Error in expr_subset: {e}")
        raise HTTPException(status_code=500, detail=str(e))
    
    # Cache as compressed bytes
    set_compressed(cache_key, arr.tobytes(), ex=3600)
    return {"data": b64_f32(arr)}

@app.post("/expr_batch")
def expr_batch(req: BatchRequest):
    """
    FAST PARQUET: Batch expression matrix.
    Lightning fast gene retrieval from Parquet!
    """
    h = _hash_idx(req.idx0)
    gene_hash = hashlib.md5(('|'.join(req.genes)).encode()).hexdigest()
    cache_key = f"expr_batch:{req.dataset_id}:{req.measurement}:{req.layer}:{h}:{gene_hash}"
    
    # Try binary cache for matrix data
    cached_bytes = get_compressed(f"{cache_key}:matrix")
    cached_meta = get_json(f"{cache_key}:meta")
    
    if cached_bytes is not None and cached_meta is not None:
        # Reconstruct from cache
        n_obs = cached_meta["n_obs"]
        n_genes = cached_meta["n_genes"]
        mat = np.frombuffer(cached_bytes, dtype=np.float32).reshape(n_obs, n_genes)
        return {
            "genes": cached_meta["genes"],
            "n_obs": n_obs,
            "n_genes": n_genes,
            "data": b64_f32(mat.reshape(-1))
        }
    
    # Get data from Parquet backend
    backend = get_backend(req.dataset_id)
    idx = np.asarray(req.idx0, dtype=np.int64) if req.idx0 else None
    
    try:
        mat = backend.get_genes_expression(req.genes, cell_indices=idx)
    except Exception:
        return {"data": None}
    
    # Cache matrix as compressed bytes and metadata as JSON
    set_compressed(f"{cache_key}:matrix", mat.tobytes(), ex=1800)
    set_json(f"{cache_key}:meta", {
        "genes": req.genes,
        "n_obs": mat.shape[0],
        "n_genes": mat.shape[1]
    }, ex=1800)
    
    return {
        "genes": req.genes,
        "n_obs": mat.shape[0],
        "n_genes": mat.shape[1],
        "data": b64_f32(mat.reshape(-1))
    }
