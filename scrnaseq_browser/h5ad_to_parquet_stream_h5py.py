#!/usr/bin/env python3
"""
Stream a large .h5ad to minimal Parquet without loading dense layers into RAM.

Supports:
- obs -> obs.parquet
- var -> var.parquet
- obsm numeric 2D -> obsm/<key>.parquet
- X (CSR) -> expression/X/chunk_XXXX.parquet (LONG: cell_id, gene, value)
- layers:
    * CSR group like raw_counts -> LONG
    * dense dataset like unlogged_normalized -> WIDE by gene chunks (cell_id + gene columns)

Usage:
python h5ad_to_parquet_stream_h5py.py \
  --h5ad BoneMarrow_May2024_outliers_removed_Palantir.h5ad \
  --outdir parquet_bm \
  --matrices X layers/raw_counts layers/unlogged_normalized layers/lognorm_pseudocount.1 \
  --obsm X_umap X_pca X_diffmap DM_EigenVectors DM_EigenVectors_multiscaled X_draw_graph_fa \
  --chunk-genes 256 \
  --compression zstd
"""

import argparse
import json
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import pandas as pd
import h5py

import pyarrow as pa
import pyarrow.parquet as pq
from scipy import sparse


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_parquet(df: pd.DataFrame, path: Path, compression: str = "zstd") -> None:
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path.as_posix(), compression=compression)


def read_h5ad_index(f: h5py.File, group: str) -> np.ndarray:
    """
    Read the index vector from /obs/_index or /var/_index or for your case:
    - obs index column: '_index'
    - var index column: 'gene'
    """
    g = f[group]
    # Common patterns:
    # - anndata stores index in "_index"
    # - or stores via "index" attr
    if "_index" in g:
        idx = g["_index"][()]
        return idx.astype(str)
    # Fallback: try attrs
    idx_name = g.attrs.get("_index", None)
    if idx_name and idx_name in g:
        return g[idx_name][()].astype(str)
    raise KeyError(f"Could not find index for group '{group}'. Keys: {list(g.keys())}")


def read_str_dataset(ds: h5py.Dataset) -> np.ndarray:
    arr = ds[()]
    # h5py may return bytes
    if arr.dtype.kind in ("S", "O"):
        return np.array([x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) for x in arr], dtype=object)
    return arr.astype(str)


def read_categorical_from_group(grp: h5py.Group) -> np.ndarray:
    """
    Read a categorical variable stored as a group with 'codes' and 'categories'.
    AnnData stores categoricals as:
      grp/categories: array of category labels
      grp/codes: integer array mapping each cell to a category index
    """
    if "categories" in grp and "codes" in grp:
        categories = grp["categories"][()]
        codes = grp["codes"][()]
        # Decode bytes if needed
        if categories.dtype.kind in ("S", "O"):
            categories = np.array([
                x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x) 
                for x in categories
            ], dtype=object)
        # Map codes to category values (-1 typically means NA/missing)
        result = np.empty(len(codes), dtype=object)
        for i, c in enumerate(codes):
            if c < 0:
                result[i] = None  # or "NA" if you prefer
            else:
                result[i] = categories[c]
        return result
    return None


def read_obs_df(f: h5py.File) -> pd.DataFrame:
    obs = f["obs"]
    cell_ids = read_str_dataset(obs["_index"]) if "_index" in obs else read_h5ad_index(f, "obs")
    data = {"cell_id": cell_ids}

    for k in obs.keys():
        if k == "_index":
            continue
        ds = obs[k]
        if isinstance(ds, h5py.Dataset):
            if ds.dtype.kind in ("S", "O"):
                data[k] = read_str_dataset(ds)
            else:
                data[k] = ds[()]
        elif isinstance(ds, h5py.Group):
            # Handle categorical variables stored as groups
            cat_values = read_categorical_from_group(ds)
            if cat_values is not None:
                data[k] = cat_values
            else:
                print(f"  ⚠️ Skipping obs/{k}: group without categories/codes")
        else:
            continue

    return pd.DataFrame(data)


def read_var_df(f: h5py.File) -> Tuple[pd.DataFrame, List[str]]:
    var = f["var"]
    # In your file, the index column is 'gene'
    if "gene" in var and isinstance(var["gene"], h5py.Dataset):
        genes = read_str_dataset(var["gene"])
    else:
        genes = read_h5ad_index(f, "var").astype(str)

    data = {"gene": genes}
    for k in var.keys():
        if k == "gene":
            continue
        ds = var[k]
        if isinstance(ds, h5py.Dataset):
            if ds.dtype.kind in ("S", "O"):
                data[k] = read_str_dataset(ds)
            else:
                data[k] = ds[()]
        else:
            continue

    return pd.DataFrame(data), list(genes.astype(str))


def export_obsm(f: h5py.File, outdir: Path, obsm_keys: List[str], cell_ids: np.ndarray, compression: str) -> List[str]:
    exported = []
    if "obsm" not in f:
        return exported
    obsm = f["obsm"]
    obsm_dir = outdir / "obsm"
    ensure_dir(obsm_dir)

    for k in obsm_keys:
        if k not in obsm:
            continue
        obj = obsm[k]
        if not isinstance(obj, h5py.Dataset):
            # groups like branch_masks/palantir_fate_probabilities -> skip
            continue
        X = obj
        if len(X.shape) != 2 or X.shape[0] != len(cell_ids):
            continue
        # numeric only
        if X.dtype.kind not in ("f", "i", "u"):
            continue

        arr = X[()].astype(np.float32, copy=False)
        cols = [f"{k}_{i}" for i in range(arr.shape[1])]
        df = pd.DataFrame(arr, columns=cols)
        df.insert(0, "cell_id", cell_ids.astype(str))
        write_parquet(df, obsm_dir / f"{k}.parquet", compression=compression)
        exported.append(k)

    return exported


def read_csr_from_group(g: h5py.Group) -> sparse.csr_matrix:
    data = g["data"][()]
    indices = g["indices"][()]
    indptr = g["indptr"][()]
    shape = tuple(g.attrs["shape"]) if "shape" in g.attrs else None
    if shape is None:
        # try infer
        # NOTE: for CSR, n_rows = len(indptr)-1, n_cols = max(indices)+1 (unsafe)
        n_rows = len(indptr) - 1
        n_cols = int(indices.max()) + 1 if indices.size > 0 else 0
        shape = (n_rows, n_cols)
    return sparse.csr_matrix((data, indices, indptr), shape=shape)


def export_sparse_long_chunks(
    csr: sparse.csr_matrix,
    mat_name: str,
    outdir: Path,
    cell_ids: np.ndarray,
    genes: List[str],
    chunk_genes: int,
    compression: str,
    value_dtype: np.dtype = np.float32,
) -> Dict[str, int]:
    expr_dir = outdir / "expression" / mat_name
    ensure_dir(expr_dir)

    gene_to_chunk: Dict[str, int] = {}
    n_genes = len(genes)

    chunk_idx = 0
    for g0 in range(0, n_genes, chunk_genes):
        g1 = min(g0 + chunk_genes, n_genes)
        chunk_gene_names = genes[g0:g1]

        sub = csr[:, g0:g1].tocoo(copy=False)
        if sub.nnz == 0:
            df = pd.DataFrame({"cell_id": [], "gene": [], "value": []})
        else:
            df = pd.DataFrame({
                "cell_id": pd.Categorical(cell_ids[sub.row].astype(str)),
                "gene": pd.Categorical([chunk_gene_names[j] for j in sub.col]),
                "value": sub.data.astype(value_dtype, copy=False),
            })
        write_parquet(df, expr_dir / f"chunk_{chunk_idx:04d}.parquet", compression=compression)

        for g in chunk_gene_names:
            gene_to_chunk[g] = chunk_idx
        chunk_idx += 1

    return gene_to_chunk


def export_dense_wide_chunks(
    ds: h5py.Dataset,
    mat_name: str,
    outdir: Path,
    cell_ids: np.ndarray,
    genes: List[str],
    chunk_genes: int,
    compression: str,
    cast_dtype: np.dtype = np.float32,
) -> Dict[str, int]:
    expr_dir = outdir / "expression" / mat_name
    ensure_dir(expr_dir)

    gene_to_chunk: Dict[str, int] = {}
    n_genes = len(genes)

    chunk_idx = 0
    for g0 in range(0, n_genes, chunk_genes):
        g1 = min(g0 + chunk_genes, n_genes)
        chunk_gene_names = genes[g0:g1]

        # KEY: slice HDF5 by columns to avoid reading whole matrix
        block = ds[:, g0:g1].astype(cast_dtype, copy=False)

        df = pd.DataFrame(block, columns=chunk_gene_names)
        df.insert(0, "cell_id", cell_ids.astype(str))
        write_parquet(df, expr_dir / f"chunk_{chunk_idx:04d}.parquet", compression=compression)

        for g in chunk_gene_names:
            gene_to_chunk[g] = chunk_idx
        chunk_idx += 1

    return gene_to_chunk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--matrices", nargs="*", default=["X", "layers/raw_counts"])
    ap.add_argument("--obsm", nargs="*", default=["X_umap", "X_pca"])
    ap.add_argument("--chunk-genes", type=int, default=256)
    ap.add_argument("--compression", default="zstd")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    ensure_dir(outdir)

    with h5py.File(args.h5ad, "r") as f:
        # obs/var
        obs_df = read_obs_df(f)
        write_parquet(obs_df, outdir / "obs.parquet", compression=args.compression)
        cell_ids = obs_df["cell_id"].to_numpy()

        var_df, genes = read_var_df(f)
        write_parquet(var_df, outdir / "var.parquet", compression=args.compression)

        exported_obsm = export_obsm(f, outdir, args.obsm, cell_ids, compression=args.compression)

        manifest = {
            "source_h5ad": Path(args.h5ad).name,
            "n_obs": int(len(cell_ids)),
            "n_vars": int(len(genes)),
            "chunk_genes": int(args.chunk_genes),
            "compression": args.compression,
            "files": {"obs": "obs.parquet", "var": "var.parquet", "obsm_dir": "obsm/", "expression_dir": "expression/"},
            "obsm_exported": exported_obsm,
            "matrices": {},
        }

        # Export matrices
        for m in args.matrices:
            if m == "X":
                Xg = f["X"]
                csr = read_csr_from_group(Xg)
                gene_to_chunk = export_sparse_long_chunks(
                    csr, "X", outdir, cell_ids, genes, args.chunk_genes, args.compression, value_dtype=np.float32
                )
                manifest["matrices"]["X"] = {
                    "format": "sparse_long_gene_chunks",
                    "chunk_files_pattern": "expression/X/chunk_{chunk:04d}.parquet",
                    "columns": ["cell_id", "gene", "value"],
                    "gene_to_chunk": gene_to_chunk,
                }
                continue

            # layers/<name>
            if m.startswith("layers/"):
                layer_name = m.split("/", 1)[1]
                obj = f["layers"][layer_name]

                if isinstance(obj, h5py.Group):
                    # CSR-like
                    csr = read_csr_from_group(obj)
                    gene_to_chunk = export_sparse_long_chunks(
                        csr, f"layers__{layer_name}", outdir, cell_ids, genes, args.chunk_genes, args.compression, value_dtype=np.float32
                    )
                    manifest["matrices"][f"layers/{layer_name}"] = {
                        "format": "sparse_long_gene_chunks",
                        "chunk_files_pattern": f"expression/layers__{layer_name}/chunk_{{chunk:04d}}.parquet",
                        "columns": ["cell_id", "gene", "value"],
                        "gene_to_chunk": gene_to_chunk,
                    }
                elif isinstance(obj, h5py.Dataset):
                    # dense dataset
                    gene_to_chunk = export_dense_wide_chunks(
                        obj, f"layers__{layer_name}", outdir, cell_ids, genes, args.chunk_genes, args.compression, cast_dtype=np.float32
                    )
                    manifest["matrices"][f"layers/{layer_name}"] = {
                        "format": "dense_wide_gene_chunks",
                        "chunk_files_pattern": f"expression/layers__{layer_name}/chunk_{{chunk:04d}}.parquet",
                        "cell_id_col": "cell_id",
                        "gene_to_chunk": gene_to_chunk,
                    }
                else:
                    raise TypeError(f"Unsupported layers/{layer_name} type: {type(obj)}")
                continue

            raise ValueError(f"Unknown matrix spec '{m}'. Use 'X' or 'layers/<name>'.")

        with open(outdir / "manifest.json", "w", encoding="utf-8") as fh:
            json.dump(manifest, fh, indent=2)

    print(f"[OK] Exported to {outdir.resolve()}")
    print(f"  - obs.parquet / var.parquet")
    print(f"  - obsm exported: {exported_obsm}")
    print(f"  - matrices: {list(manifest['matrices'].keys())}")


if __name__ == "__main__":
    main()
