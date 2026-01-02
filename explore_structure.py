import h5py
import os
import sys

# File to inspect
file_path = 'BoneMarrow_May2024_outliers_removed_Palantir.h5ad'

def print_structure(name, obj):
    """
    Callback function for h5py.visititems to print structure.
    """
    indent = "  " * (name.count('/') + 1)
    base_name = name.split('/')[-1]
    
    if isinstance(obj, h5py.Dataset):
        print(f"{indent}- {base_name} (Dataset): {obj.shape} - {obj.dtype}")
    elif isinstance(obj, h5py.Group):
        print(f"{indent}+ {base_name} (Group)")
        # Special handling for obs and var to list columns without recursing too deep if it's just columns
        if base_name in ['obs', 'var']:
            print(f"{indent}  [Columns]: {list(obj.keys())}")

def explore_h5ad_structure():
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    print(f"Inspecting structure of '{file_path}' using h5py...")
    print("=" * 60)

    try:
        with h5py.File(file_path, 'r') as f:
            # Print top-level keys
            print(f"Top-level keys: {list(f.keys())}")
            print("-" * 60)
            
            # Inspect specific important groups
            for key in ['obs', 'var', 'obsm', 'layers', 'uns']:
                if key in f:
                    print(f"\n[{key}]")
                    group = f[key]
                    if isinstance(group, h5py.Group):
                        keys = list(group.keys())
                        print(f"  Count: {len(keys)}")
                        print(f"  Keys: {keys[:20]} {'...' if len(keys) > 20 else ''}")
                        
                        # If it's obs or var, we might want to see the index
                        if '_index' in group.attrs:
                            print(f"  Index column: {group.attrs['_index']}")
                    elif isinstance(group, h5py.Dataset):
                         print(f"  Dataset: {group.shape}")
                else:
                    print(f"\n[{key}] - Not found")

            print("\n" + "=" * 60)
            print("Detailed Structure (First 2 levels):")
            
            def visit_func(name, node):
                # Limit depth to avoid spamming
                if name.count('/') < 2:
                    indent = "  " * name.count('/')
                    if isinstance(node, h5py.Group):
                        print(f"{indent}+ {name.split('/')[-1]}/")
                    else:
                        print(f"{indent}- {name.split('/')[-1]} {node.shape} {node.dtype}")

            f.visititems(visit_func)

    except Exception as e:
        print(f"Error reading file: {e}")

if __name__ == "__main__":
    explore_h5ad_structure()
