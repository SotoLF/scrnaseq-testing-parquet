#!/usr/bin/env python3
"""
Add missing categorical columns from h5ad to existing obs.parquet.
Uses h5py to read only the categorical variables without loading the entire object.

Usage:
    python add_categorical_to_obs.py \
        --h5ad BoneMarrow_May2024_outliers_removed_Palantir.h5ad \
        --parquet-dir parquet_bm
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import h5py
import pyarrow as pa
import pyarrow.parquet as pq


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
                result[i] = None
            else:
                result[i] = categories[c]
        return result
    return None


def get_categorical_columns(f: h5py.File) -> dict:
    """Extract only categorical columns from obs group."""
    obs = f["obs"]
    categoricals = {}
    
    for k in obs.keys():
        if k == "_index":
            continue
        item = obs[k]
        if isinstance(item, h5py.Group):
            # This is likely a categorical
            cat_values = read_categorical_from_group(item)
            if cat_values is not None:
                categoricals[k] = cat_values
                print(f"  ✓ Found categorical: {k} ({len(np.unique(cat_values))} unique values)")
    
    return categoricals


def main():
    parser = argparse.ArgumentParser(description="Add categorical columns to obs.parquet")
    parser.add_argument("--h5ad", required=True, help="Path to h5ad file")
    parser.add_argument("--parquet-dir", required=True, help="Path to parquet output directory")
    parser.add_argument("--backup", action="store_true", help="Create backup of original obs.parquet")
    args = parser.parse_args()
    
    h5ad_path = Path(args.h5ad)
    parquet_dir = Path(args.parquet_dir)
    obs_path = parquet_dir / "obs.parquet"
    
    if not h5ad_path.exists():
        print(f"❌ h5ad file not found: {h5ad_path}")
        return
    
    if not obs_path.exists():
        print(f"❌ obs.parquet not found: {obs_path}")
        return
    
    print(f"📖 Reading categorical columns from {h5ad_path}...")
    with h5py.File(h5ad_path, "r") as f:
        categoricals = get_categorical_columns(f)
    
    if not categoricals:
        print("⚠️ No categorical columns found in h5ad")
        return
    
    print(f"\n📖 Reading existing obs.parquet...")
    obs_df = pd.read_parquet(obs_path)
    print(f"   Current columns: {list(obs_df.columns)}")
    print(f"   Shape: {obs_df.shape}")
    
    # Check which categoricals are already present
    new_cols = []
    for col_name, col_data in categoricals.items():
        if col_name in obs_df.columns:
            print(f"  ⚠️ Column '{col_name}' already exists, skipping")
        else:
            obs_df[col_name] = col_data
            new_cols.append(col_name)
            print(f"  ✓ Added column: {col_name}")
    
    if not new_cols:
        print("\n✅ No new columns to add, obs.parquet is already complete")
        return
    
    # Backup original if requested
    if args.backup:
        backup_path = parquet_dir / "obs.parquet.backup"
        import shutil
        shutil.copy(obs_path, backup_path)
        print(f"\n💾 Backup saved to: {backup_path}")
    
    # Save updated obs.parquet
    print(f"\n💾 Saving updated obs.parquet...")
    table = pa.Table.from_pandas(obs_df, preserve_index=False)
    pq.write_table(table, obs_path, compression="zstd")
    
    print(f"\n✅ Done! Added {len(new_cols)} categorical columns: {new_cols}")
    print(f"   New shape: {obs_df.shape}")
    print(f"   All columns: {list(obs_df.columns)}")


if __name__ == "__main__":
    main()
