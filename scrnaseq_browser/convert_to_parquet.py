#!/usr/bin/env python3
"""
Export minimal AnnData for a Shiny app:
- obs metadata (cells)  -> parquet
- var metadata (genes)  -> parquet
- selected embeddings in obsm -> parquet
- selected count matrices (layers and/or X) -> parquet, chunked by genes (wide format)

Why wide-by-gene chunks?
- Parquet is columnar: you can read only one gene column efficiently.
- Chunking avoids a single huge file and keeps reads bounded.

Example:
python adata_to_parquet_min.py \
  --h5ad input.h5ad \
  --outdir parquet_out \
  --layers raw cellbender MAGIC_imputed_data \
  --obsm X_umap X_pca \
  --chunk-genes 512 \
  --compression zstd
"""

import argparse
import json
import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

try:
    import anndata as ad
except ImportError as e:
    raise SystemExit("Missing dependency: anndata. Install with: pip install anndata") from e

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError as e:
    raise SystemExit("Missing dependency: pyarrow. Install with: pip install pyarrow") from e

from scipy import sparse


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def write_parquet(df: pd.DataFrame, path: Path, compression: str = "zstd") -> None:
    table = pa.Table.from_pandas(df, preserve_index=False)
    pq.write_table(table, path.as_posix(), compression=compression)


def get_matrix(adata, layer: str):
    """
    layer can be:
      - "X" to export adata.X
      - a layer name in adata.layers
    Returns a matrix-like object with shape (n_obs, n_vars).
    """
    if layer == "X":
        return adata.X
    if layer in adata.layers:
        return adata.layers[layer]
    raise ValueError(f"Layer '{layer}' not found (use 'X' or one of adata.layers keys).")


def slice_to_dense(mat, rows: slice, cols: slice, dtype=np.float32) -> np.ndarray:
    """
    Convert a (rows, cols) slice to a dense numpy array.
    Handles numpy arrays, scipy sparse matrices, and array-like.
    """
    sub = mat[rows, cols]
    if sparse.issparse(sub):
        arr = sub.toarray()
    else:
        arr = np.asarray(sub)
    if dtype is not None:
        arr = arr.astype(dtype, copy=False)
    return arr


def export_obs_var(adata, outdir: Path, compression: str) -> None:
    # obs
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", obs.index.astype(str))
    write_parquet(obs.reset_index(drop=True), outdir / "obs.parquet", compression=compression)

    # var
    var = adata.var.copy()
    var.insert(0, "gene", var.index.astype(str))
    write_parquet(var.reset_index(drop=True), outdir / "var.parquet", compression=compression)


def export_obsm(adata, outdir: Path, obsm_keys: List[str], compression: str) -> None:
    obsm_dir = outdir / "obsm"
    ensure_dir(obsm_dir)

    cell_ids = adata.obs_names.astype(str)

    for k in obsm_keys:
        if k not in adata.obsm:
            raise ValueError(f"obsm key '{k}' not found. Available: {list(adata.obsm.keys())}")
        X = np.asarray(adata.obsm[k])
        if X.ndim != 2 or X.shape[0] != adata.n_obs:
            raise ValueError(f"obsm['{k}'] has unexpected shape {X.shape}; expected (n_obs, d).")

        cols = [f"{k}_{i}" for i in range(X.shape[1])]
        df = pd.DataFrame(X, columns=cols)
        df.insert(0, "cell_id", cell_ids)
        write_parquet(df, obsm_dir / f"{k}.parquet", compression=compression)


def export_layer_chunked_by_genes(
    adata,
    layer: str,
    outdir: Path,
    chunk_genes: int,
    compression: str,
    dtype: str,
) -> Dict[str, int]:
    """
    Exports one matrix (layer or X) into multiple parquet files, each containing:
      - cell_id
      - gene columns for a contiguous gene block [g0:g1)

    Returns a dict mapping gene -> chunk_index for manifest creation.
    """
    mat = get_matrix(adata, layer=layer)

    # Decide dtype
    # Counts layers often should be int32; imputed/normalized often float32.
    np_dtype = np.int32 if dtype == "int32" else np.float32

    layer_dir = outdir / "X_layers" / layer
    ensure_dir(layer_dir)

    cell_ids = adata.obs_names.astype(str)
    genes = adata.var_names.astype(str)

    n_cells, n_genes = adata.n_obs, adata.n_vars
    gene_to_chunk: Dict[str, int] = {}

    # Export chunks
    chunk_idx = 0
    for g0 in range(0, n_genes, chunk_genes):
        g1 = min(g0 + chunk_genes, n_genes)

        # Build dense block (n_cells x (g1-g0))
        block = slice_to_dense(mat, slice(None), slice(g0, g1), dtype=np_dtype)

        # DataFrame with gene columns
        chunk_genes_names = genes[g0:g1].tolist()
        df = pd.DataFrame(block, columns=chunk_genes_names)
        df.insert(0, "cell_id", cell_ids)

        # Write
        chunk_path = layer_dir / f"chunk_{chunk_idx:03d}.parquet"
        write_parquet(df, chunk_path, compression=compression)

        # Update manifest mapping
        for g in chunk_genes_names:
            gene_to_chunk[g] = chunk_idx

        chunk_idx += 1

    return gene_to_chunk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input AnnData .h5ad path")
    ap.add_argument("--outdir", required=True, help="Output directory for parquet dataset")
    ap.add_argument(
        "--layers",
        nargs="*",
        default=["raw", "cellbender", "MAGIC_imputed_data"],
        help="Which matrices to export. Use 'X' for adata.X, or any key in adata.layers",
    )
    ap.add_argument(
        "--obsm",
        nargs="*",
        default=["X_umap", "X_pca"],
        help="Which embeddings to export from adata.obsm",
    )
    ap.add_argument("--chunk-genes", type=int, default=512, help="Genes per parquet chunk file")
    ap.add_argument("--compression", default="zstd", help="Parquet compression: zstd, snappy, gzip, none")
    ap.add_argument(
        "--dtype",
        default="float32",
        choices=["float32", "int32"],
        help="Storage dtype for expression matrices",
    )
    args = ap.parse_args()

    h5ad = Path(args.h5ad)
    outdir = Path(args.outdir)
    ensure_dir(outdir)

    # Read
    adata = ad.read_h5ad(h5ad.as_posix())

    # Export minimal pieces
    export_obs_var(adata, outdir, compression=args.compression)
    export_obsm(adata, outdir, args.obsm, compression=args.compression)

    # Export layers (chunked)
    manifest = {
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "cell_id_col": "cell_id",
        "genes": adata.var_names.astype(str).tolist(),
        "cells": adata.obs_names.astype(str).tolist(),
        "chunk_genes": int(args.chunk_genes),
        "compression": args.compression,
        "layers": {},
        "obsm": args.obsm,
        "files": {
            "obs": "obs.parquet",
            "var": "var.parquet",
            "obsm_dir": "obsm/",
            "layers_dir": "X_layers/",
        },
    }

    for layer in args.layers:
        gene_to_chunk = export_layer_chunked_by_genes(
            adata=adata,
            layer=layer,
            outdir=outdir,
            chunk_genes=args.chunk_genes,
            compression=args.compression,
            dtype=args.dtype,
        )
        manifest["layers"][layer] = {
            "format": "parquet_wide_gene_chunks",
            "dtype": args.dtype,
            "gene_to_chunk": gene_to_chunk,
            "chunk_files_pattern": f"X_layers/{layer}/chunk_{{chunk:03d}}.parquet",
        }

    # Write manifest
    with open(outdir / "manifest.json", "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"[OK] Exported to: {outdir.resolve()}")
    print("Files:")
    print(f"  - {outdir/'obs.parquet'}")
    print(f"  - {outdir/'var.parquet'}")
    print(f"  - {outdir/'obsm/'}")
    print(f"  - {outdir/'X_layers/'}")
    print(f"  - {outdir/'manifest.json'}")


if __name__ == "__main__":
    main()