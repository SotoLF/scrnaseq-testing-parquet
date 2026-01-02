import scanpy as sc
import os

# Define the file path
file_path = 'strict_epdsc_annotated_data.h5'
file_path = 'BoneMarrow_May2024_outliers_removed_Palantir.h5ad'

def inspect_h5_file():
    # Check if the file exists in the current directory
    if not os.path.exists(file_path):
        print(f"Error: The file '{file_path}' was not found in the current directory.")
        return

    print(f"Attempting to load '{file_path}' with scanpy...")

    try:
        # Read the file
        if file_path.endswith('.h5ad'):
            try:
                print("Attempting to load fully into memory...")
                adata = sc.read_h5ad(file_path)
            except Exception as e:
                print(f"\n⚠️ Full load failed: {e}")
                print("Attempting to load in 'backed' mode (memory efficient)...")
                adata = sc.read_h5ad(file_path, backed='r')
        else:
            adata = sc.read(file_path)

        print("\n" + "="*40)
        print("       AnnData Object Summary")
        print("="*40)
        print(adata)
        
        print("\n" + "-"*40)
        print("1. Observations (Cells) - adata.obs")
        print("-"*40)
        print(f"Shape: {adata.obs.shape}")
        print("First 5 rows:")
        print(adata.obs.head())

        print("\n" + "-"*40)
        print("2. Variables (Genes) - adata.var")
        print("-"*40)
        print(f"Shape: {adata.var.shape}")
        print("First 5 rows:")
        print(adata.var.head())

        if hasattr(adata, 'obsm') and adata.obsm:
            print("\n" + "-"*40)
            print("3. Embeddings/Matrices - adata.obsm")
            print("-"*40)
            print(f"Keys: {list(adata.obsm.keys())}")

        if hasattr(adata, 'uns') and adata.uns:
            print("\n" + "-"*40)
            print("4. Unstructured Data - adata.uns")
            print("-"*40)
            print(f"Keys: {list(adata.uns.keys())}")

        if hasattr(adata, 'layers') and adata.layers:
            print("\n" + "-"*40)
            print("5. Layers - adata.layers")
            print("-"*40)
            print(f"Keys: {list(adata.layers.keys())}")

    except Exception as e:
        print(f"\nError loading file: {e}")
        print("Tip: If this is a raw 10x HDF5 file, you might need to use 'sc.read_10x_h5'.")

if __name__ == "__main__":
    inspect_h5_file()
