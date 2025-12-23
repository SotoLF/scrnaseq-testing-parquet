"""Fast Parquet backend for scRNA-seq data.

Design goals:
- Keep startup fast.
- Keep RAM usage low (avoid loading whole matrices/metadata eagerly).
- Read only the needed Parquet columns per request.
"""

import json
import os
from pathlib import Path
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


_CHUNK_CACHE_SIZE = int(os.getenv("CHUNK_CACHE_SIZE", "2"))


class ParquetBackend:
    """
    Optimized Parquet reader for gene expression data.
    Provides instant gene retrieval (~10-50ms vs 5-20s from SOMA).
    """
    
    def __init__(self, parquet_dir: str = None):
        if parquet_dir is None:
            # Check environment variable or use S3 path
            parquet_dir = os.getenv('PARQUET_DATA_DIR', 's3://luis-scrnaseq-data/parquet_data')
        
        # Convert relative path to absolute if needed
        if not parquet_dir.startswith('s3://') and not os.path.isabs(parquet_dir):
            # Try to use PROJECT_ROOT env var first, otherwise calculate from file location
            project_root = os.getenv('PROJECT_ROOT')
            if project_root:
                parquet_dir = os.path.join(project_root, parquet_dir.lstrip('./'))
            else:
                # Get project root (2 levels up from backend/app/)
                # Use resolve() to get absolute path
                backend_app_dir = Path(__file__).resolve().parent
                project_root = backend_app_dir.parent.parent
                parquet_dir = str(project_root / parquet_dir.lstrip('./'))
            
            # Ensure the path exists
            if not os.path.exists(parquet_dir):
                raise FileNotFoundError(
                    f"Parquet data directory not found: {parquet_dir}\n"
                    f"Make sure PARQUET_DATA_DIR in .env points to the correct location."
                )
        
        self.parquet_dir = parquet_dir
        self._load_metadata()
        
        print(f"✅ Initialized ParquetBackend from {parquet_dir}")
        print(f"   - {self.n_cells:,} cells")
        print(f"   - {self.n_genes:,} genes")
        print(f"   - {len(self.available_layers)} layers: {self.available_layers}")
    
    def _load_metadata(self):
        """Load manifest and metadata tables"""
        # Load manifest
        manifest_path = f"{self.parquet_dir}/manifest.json"
        with pd.io.common.get_handle(manifest_path, 'r') as handle:
            self.manifest = json.load(handle.handle)
        
        # Support new manifest format from convert_to_parquet.py
        if 'n_obs' in self.manifest:
            self.n_cells = self.manifest['n_obs']
            self.n_genes = self.manifest['n_vars']
            
            # Handle layers/matrices
            if 'layers' in self.manifest:
                self.layers_info = self.manifest['layers']
            elif 'matrices' in self.manifest:
                self.layers_info = self.manifest['matrices']
            else:
                self.layers_info = {}
            
            self.available_layers = list(self.layers_info.keys())
            
            # Build layer alias map for flexible lookup
            # Maps common names like "X" to actual layer keys
            self._layer_alias = {}
            for layer_key in self.available_layers:
                # Map exact name
                self._layer_alias[layer_key] = layer_key
                # Map simplified names (e.g., "raw_counts" for "layers/raw_counts")
                if '/' in layer_key:
                    simple_name = layer_key.split('/')[-1]
                    if simple_name not in self._layer_alias:
                        self._layer_alias[simple_name] = layer_key
            
            # If no "X" layer exists, map "X" to the first available layer
            if "X" not in self._layer_alias and self.available_layers:
                self._layer_alias["X"] = self.available_layers[0]
                print(f"   Note: No 'X' layer found, mapping 'X' to '{self.available_layers[0]}'")
            
            # Handle obsm
            if 'obsm' in self.manifest:
                self.available_embeddings = self.manifest['obsm']
            elif 'obsm_exported' in self.manifest:
                self.available_embeddings = self.manifest['obsm_exported']
            else:
                self.available_embeddings = []
            
            # Use the first layer to establish gene -> chunk mapping
            if self.available_layers:
                first_layer = self.available_layers[0]
                self.gene_to_chunk = self.layers_info[first_layer].get('gene_to_chunk', {})
            else:
                self.gene_to_chunk = {}
        else:
            # Fallback to old format
            self.n_cells = self.manifest['n_cells']
            self.n_genes = self.manifest['n_genes']
            self.available_layers = ["X"]
            self.available_embeddings = ["X_umap", "X_pca"] # Default guess
            
            # Build gene -> chunk mapping
            self.gene_to_chunk = {}
            for chunk in self.manifest['chunks']:
                for gene in chunk['genes']:
                    self.gene_to_chunk[gene] = chunk['chunk_id']
        
        # Paths used for Parquet reads.
        self._genes_path = f"{self.parquet_dir}/var.parquet"
        self._cells_path = f"{self.parquet_dir}/obs.parquet"
        
        # Load cell IDs for mapping (needed for sparse_long format)
        self.cell_ids = None
        self.cell_id_to_idx = None
        try:
            # Try to read cell_id column
            cell_id_col = self.manifest.get("cell_id_col", "cell_id")
            try:
                tbl = pq.read_table(self._cells_path, columns=[cell_id_col])
                self.cell_ids = tbl.column(0).to_numpy().astype(str)
            except Exception:
                # Fallback to index if column not found
                tbl = pq.read_table(self._cells_path)
                self.cell_ids = tbl.column(0).to_numpy().astype(str)
            
            self.cell_id_to_idx = {cid: i for i, cid in enumerate(self.cell_ids)}
        except Exception as e:
            print(f"Warning: Could not load cell IDs: {e}")

        # Handle obsm dir from manifest if available
        if 'files' in self.manifest and 'obsm_dir' in self.manifest['files']:
             # Strip trailing slash from obsm_dir to avoid double slashes later
             obsm_rel = self.manifest['files']['obsm_dir'].rstrip('/')
             self._obsm_dir = f"{self.parquet_dir}/{obsm_rel}"
        else:
             self._obsm_dir = f"{self.parquet_dir}/obsm"

        # Preload metadata into RAM (fast) only for local paths.
        preload = os.getenv("PRELOAD_METADATA", "true").lower() == "true"
        is_s3 = str(self.parquet_dir).startswith("s3://")

        self.gene_metadata: Optional[pd.DataFrame] = None
        self.cell_metadata: Optional[pd.DataFrame] = None
        self.embeddings_cache: Dict[str, pd.DataFrame] = {}

        gene_ids = None
        try:
            if preload and not is_s3:
                self.gene_metadata = pd.read_parquet(self._genes_path)
                if "gene" in self.gene_metadata.columns:
                    gene_ids = self.gene_metadata["gene"].astype(str).tolist()
                elif "gene_id" in self.gene_metadata.columns:
                    gene_ids = self.gene_metadata["gene_id"].astype(str).tolist()
            else:
                # Try reading 'gene' column first
                try:
                    gene_tbl = pq.read_table(self._genes_path, columns=["gene"])
                except:
                    gene_tbl = pq.read_table(self._genes_path, columns=["gene_id"])
                gene_ids = [str(x) for x in gene_tbl.column(0).to_pylist()]
        except Exception:
            gene_ids = None

        try:
            if preload and not is_s3:
                self.cell_metadata = pd.read_parquet(self._cells_path)
        except Exception:
            self.cell_metadata = None

        # Map internal Parquet gene columns (gene_0...) to display names (prefer gene_id).
        self.internal_to_display: Dict[str, str] = {}
        self.display_to_internal: Dict[str, str] = {}
        if gene_ids is not None:
            for i, display in enumerate([str(x) for x in gene_ids]):
                # In new format, column name IS the gene name.
                # So internal == display.
                self.internal_to_display[display] = display
                self.display_to_internal[display] = display
    
    def _resolve_layer(self, layer: str) -> str:
        """Resolve a layer name to the actual key in layers_info.
        
        This handles cases where:
        - User requests 'X' but we have 'raw' as first layer
        - User requests 'raw_counts' but we have 'layers/raw_counts'
        """
        if hasattr(self, '_layer_alias') and layer in self._layer_alias:
            return self._layer_alias[layer]
        # Fallback: return as-is
        return layer
    
    def _chunk_path(self, chunk_id: int, layer: str = "X") -> str:
        # Resolve layer name to actual key
        resolved_layer = self._resolve_layer(layer)
        
        # Check if we have a specific pattern for this layer in manifest
        if hasattr(self, 'layers_info') and resolved_layer in self.layers_info:
            pattern = self.layers_info[resolved_layer].get('chunk_files_pattern')
            if pattern:
                # pattern is like "X_layers/raw/chunk_{chunk:03d}.parquet"
                # or "expression/X/chunk_{chunk:04d}.parquet"
                rel_path = pattern.format(chunk=chunk_id)
                return f"{self.parquet_dir}/{rel_path}"
        
        # Fallback for old format or missing pattern
        return f"{self.parquet_dir}/X_layers/{layer}/chunk_{chunk_id:03d}.parquet"

    @lru_cache(maxsize=_CHUNK_CACHE_SIZE)
    def _load_chunk_df(self, chunk_id: int, layer: str = "X") -> pd.DataFrame:
        """Load a full chunk into RAM (small LRU cache for speed)."""
        return pd.read_parquet(self._chunk_path(chunk_id, layer))

    def _read_chunk_columns(self, chunk_id: int, columns: List[str], layer: str = "X") -> pa.Table:
        """Read only requested columns from a chunk (low RAM, slower)."""
        return pq.read_table(self._chunk_path(chunk_id, layer), columns=columns)

    
    def get_gene_expression(self, gene_id: str, cell_indices: Optional[np.ndarray] = None, layer: str = "X") -> np.ndarray:
        """
        Get expression values for a gene.
        """
        internal_gene = gene_id 

        # Find which chunk contains this gene
        if internal_gene not in self.gene_to_chunk:
            raise ValueError(f"Gene {gene_id} not found in dataset")

        chunk_id = self.gene_to_chunk[internal_gene]
        
        # Resolve layer name to actual key
        resolved_layer = self._resolve_layer(layer)
        
        # Check format
        layer_info = self.layers_info.get(resolved_layer, {})
        fmt = layer_info.get("format", "parquet_wide_gene_chunks")

        if fmt == "sparse_long_gene_chunks":
            return self._get_gene_expression_long(chunk_id, internal_gene, cell_indices, resolved_layer)

        # Prefer cached chunk for speed; fall back to column read.
        expr = None
        if _CHUNK_CACHE_SIZE > 0 and not str(self.parquet_dir).startswith("s3://"):
            try:
                chunk_df = self._load_chunk_df(chunk_id, resolved_layer)
                expr = chunk_df[internal_gene].to_numpy(dtype=np.float32, copy=False)
            except Exception:
                expr = None
        if expr is None:
            tbl = self._read_chunk_columns(chunk_id, [internal_gene], resolved_layer)
            expr = tbl.column(0).to_numpy(zero_copy_only=False).astype(np.float32, copy=False)
        
        # Subset if requested
        if cell_indices is not None:
            expr = expr[cell_indices]
        
        return expr

    def _get_gene_expression_long(self, chunk_id: int, gene: str, cell_indices: Optional[np.ndarray], layer: str) -> np.ndarray:
        """Handle sparse long format: cell_id, gene, value"""
        # We must read the whole chunk (or filter pushdown)
        # Since chunks are small (256 genes), reading whole chunk is OK.
        # But we can try to filter by gene if pyarrow supports it efficiently.
        
        path = self._chunk_path(chunk_id, layer)
        
        # Read table with filter
        # Filters: [('gene', '==', gene)]
        try:
            tbl = pq.read_table(path, filters=[('gene', '==', gene)])
        except Exception:
            # Fallback if filters not supported or error
            tbl = pq.read_table(path)
            # Filter in pandas
        
        df = tbl.to_pandas()
        if 'gene' in df.columns:
            df = df[df['gene'] == gene]
        
        # Now we have cell_id and value.
        # We need to map cell_id to index.
        if self.cell_id_to_idx is None:
            raise ValueError("Cell IDs not loaded, cannot map sparse data")
            
        # Create result array
        n_cells = self.n_cells
        res = np.zeros(n_cells, dtype=np.float32)
        
        if df.empty:
            if cell_indices is not None:
                return res[cell_indices]
            return res
            
        # Map IDs to indices
        # df['cell_id'] might be categorical or string
        cids = df['cell_id'].astype(str).values
        vals = df['value'].values
        
        # Vectorized mapping is hard without a perfect hash or integer map.
        # But we can use the dict.
        # For speed, maybe we can rely on the fact that cell_ids are usually consistent?
        # No, sparse format skips zeros.
        
        # Slow loop?
        # indices = [self.cell_id_to_idx[cid] for cid in cids if cid in self.cell_id_to_idx]
        # This is slow in Python.
        
        # Faster: use pandas map
        # But we need a Series.
        # s = pd.Series(cids)
        # idxs = s.map(self.cell_id_to_idx).values
        
        # Even faster: if we have self.cell_ids array, we can use searchsorted if sorted?
        # No, cell_ids are not guaranteed sorted.
        
        # Let's use the dict map, it's reasonably fast for 10k-100k cells.
        # But wait, if we have 100k cells, this loop is 100k ops.
        
        # Optimization:
        # If we pre-computed a map from "cell_id string" to "int index" in a way that pandas can use?
        # pd.Index(self.cell_ids).get_indexer(cids)
        
        # This is fast!
        if not hasattr(self, '_cell_id_index'):
            self._cell_id_index = pd.Index(self.cell_ids)
            
        idxs = self._cell_id_index.get_indexer(cids)
        
        # Filter out -1 (not found)
        mask = idxs != -1
        idxs = idxs[mask]
        vals = vals[mask]
        
        res[idxs] = vals
        
        if cell_indices is not None:
            return res[cell_indices]
        return res

    
    def get_genes_expression(self, gene_ids: List[str], cell_indices: Optional[np.ndarray] = None, layer: str = "X") -> np.ndarray:
        """
        Get expression matrix for multiple genes.
        
        Args:
            gene_ids: List of gene identifiers
            cell_indices: Optional array of cell indices to subset
            layer: Layer name (default: "X")
        
        Returns:
            numpy array of shape (n_cells, n_genes)
        """
        n_cells = len(cell_indices) if cell_indices is not None else self.n_cells
        n_genes = len(gene_ids)

        # Resolve layer name to actual key
        resolved_layer = self._resolve_layer(layer)

        # Resolve display names -> internal columns.
        internal_genes: List[Optional[str]] = [g for g in gene_ids]

        # Group genes by chunk id to minimize I/O.
        by_chunk: Dict[int, List[Tuple[int, str]]] = {}
        for out_idx, internal_gene in enumerate(internal_genes):
            if internal_gene not in self.gene_to_chunk:
                continue
            chunk_id = self.gene_to_chunk[internal_gene]
            by_chunk.setdefault(chunk_id, []).append((out_idx, internal_gene))

        matrix = np.zeros((n_cells, n_genes), dtype=np.float32)
        
        # Check format using resolved layer
        layer_info = self.layers_info.get(resolved_layer, {})
        fmt = layer_info.get("format", "parquet_wide_gene_chunks")

        for chunk_id, genes_in_chunk in by_chunk.items():
            if fmt == "sparse_long_gene_chunks":
                # Handle long format
                # We can read the chunk once and filter for all genes in this chunk
                target_genes = set(g for _, g in genes_in_chunk)
                path = self._chunk_path(chunk_id, resolved_layer)
                
                # Filter pushdown for multiple genes? 
                # filters=[('gene', 'in', target_genes)] is supported in newer pyarrow
                try:
                    tbl = pq.read_table(path, filters=[('gene', 'in', target_genes)])
                except Exception:
                    tbl = pq.read_table(path)
                
                df = tbl.to_pandas()
                if 'gene' in df.columns:
                    df = df[df['gene'].isin(target_genes)]
                
                if df.empty:
                    continue
                    
                # Map IDs
                if not hasattr(self, '_cell_id_index'):
                    self._cell_id_index = pd.Index(self.cell_ids)
                
                cids = df['cell_id'].astype(str).values
                idxs = self._cell_id_index.get_indexer(cids)
                
                # We have (cell_idx, gene, value)
                # We need to fill matrix[cell_idx, out_idx]
                # Map gene to out_idx
                gene_to_out_idx = {g: i for i, g in genes_in_chunk}
                
                # Iterate? Or vectorized?
                # df['out_idx'] = df['gene'].map(gene_to_out_idx)
                # But map is slow.
                
                # Loop over genes in chunk
                for out_idx, g in genes_in_chunk:
                    sub = df[df['gene'] == g]
                    if sub.empty:
                        continue
                    
                    # Get indices for this gene
                    # We already computed idxs for the whole df, but we need to subset
                    # It's easier to do it per gene
                    
                    sub_cids = sub['cell_id'].astype(str).values
                    sub_idxs = self._cell_id_index.get_indexer(sub_cids)
                    sub_vals = sub['value'].values
                    
                    mask = sub_idxs != -1
                    valid_idxs = sub_idxs[mask]
                    valid_vals = sub_vals[mask]
                    
                    # If cell_indices provided, we need to map global idx -> local idx
                    # This is getting complicated.
                    # Simpler: fill global matrix, then subset at the end?
                    # Or: only fill if valid_idxs in cell_indices?
                    
                    if cell_indices is not None:
                        # We need to map global index to result index (0..n_cells-1)
                        # This is slow.
                        # Better: fill a full-size array then subset?
                        # If n_cells is huge (100k), full array is 400KB per gene. OK.
                        
                        full_col = np.zeros(self.n_cells, dtype=np.float32)
                        full_col[valid_idxs] = valid_vals
                        matrix[:, out_idx] = full_col[cell_indices]
                    else:
                        matrix[valid_idxs, out_idx] = valid_vals
                
                continue

            # Prefer cached full chunk for speed.
            chunk_df = None
            if _CHUNK_CACHE_SIZE > 0 and not str(self.parquet_dir).startswith("s3://"):
                try:
                    chunk_df = self._load_chunk_df(chunk_id, resolved_layer)
                except Exception:
                    chunk_df = None

            if chunk_df is not None:
                for out_idx, internal_gene in genes_in_chunk:
                    col = chunk_df[internal_gene].to_numpy(dtype=np.float32, copy=False)
                    if cell_indices is not None:
                        col = col[cell_indices]
                    matrix[:, out_idx] = col
            else:
                cols = [g for _, g in genes_in_chunk]
                tbl = self._read_chunk_columns(chunk_id, cols, resolved_layer)
                for j, (out_idx, _) in enumerate(genes_in_chunk):
                    col = tbl.column(j).to_numpy(zero_copy_only=False).astype(np.float32, copy=False)
                    if cell_indices is not None:
                        col = col[cell_indices]
                    matrix[:, out_idx] = col

        return matrix
    
    def get_obs_column(self, column: str) -> np.ndarray:
        """Get a column from cell metadata (prefer in-memory for speed)."""
        if self.cell_metadata is not None and column in self.cell_metadata.columns:
            return self.cell_metadata[column].values
        tbl = pq.read_table(self._cells_path, columns=[column])
        series = tbl.column(0).to_pandas(self_destruct=True)
        return series.values
    
    def get_obs_types(self) -> Dict[str, str]:
        """Return mapping of obs columns to type: 'num' or 'cat'."""
        out: Dict[str, str] = {}

        if self.cell_metadata is not None:
            for col in self.cell_metadata.columns:
                s = self.cell_metadata[col]
                if pd.api.types.is_numeric_dtype(s.dtype) or pd.api.types.is_bool_dtype(s.dtype):
                    out[col] = "num"
                else:
                    out[col] = "cat"
            return out

        pf = pq.ParquetFile(self._cells_path)
        schema = pf.schema_arrow
        for field in schema:
            t = field.type
            if pa.types.is_boolean(t) or pa.types.is_integer(t) or pa.types.is_floating(t):
                out[field.name] = "num"
            else:
                out[field.name] = "cat"
        return out
    
    def get_obsm_xy(self, embedding_name: str = "X_umap") -> Tuple[np.ndarray, np.ndarray]:
        """
        Get X,Y coordinates from embedding.
        
        Args:
            embedding_name: Name of embedding (e.g., 'X_umap', 'X_pca')
        
        Returns:
            (x, y) tuple of numpy arrays
        """
        # Check cache
        if embedding_name in self.embeddings_cache:
            df = self.embeddings_cache[embedding_name]
            return df.iloc[:, 0].values, df.iloc[:, 1].values

        # Read from file
        # File path: obsm/{embedding_name}.parquet
        # Ensure no double slashes
        path = f"{self._obsm_dir}/{embedding_name}.parquet"
        if path.startswith("s3://"):
             # Clean up double slashes in S3 path (except the protocol part)
             parts = path.split("s3://", 1)
             path = "s3://" + parts[1].replace("//", "/")
        else:
             path = os.path.normpath(path)
        
        try:
            # We expect columns like X_umap_0, X_umap_1
            col_x = f"{embedding_name}_0"
            col_y = f"{embedding_name}_1"
            
            # If using S3, pyarrow handles s3:// paths directly
            try:
                tbl = pq.read_table(path, columns=[col_x, col_y])
                x = tbl.column(0).to_numpy(zero_copy_only=False)
                y = tbl.column(1).to_numpy(zero_copy_only=False)
            except Exception:
                # Fallback: read all columns and find the numeric ones
                tbl = pq.read_table(path)
                names = tbl.column_names
                
                # Filter out "cell_id" or index columns
                numeric_cols = []
                for name in names:
                    if name == "cell_id" or name.startswith("__index"):
                        continue
                    # Check type
                    if pa.types.is_floating(tbl.schema.field(name).type):
                        numeric_cols.append(name)
                
                if len(numeric_cols) >= 2:
                    # Sort to ensure X, Y order (assuming suffix _0, _1 or similar)
                    numeric_cols.sort()
                    x = tbl.column(numeric_cols[0]).to_numpy(zero_copy_only=False)
                    y = tbl.column(numeric_cols[1]).to_numpy(zero_copy_only=False)
                elif tbl.num_columns >= 2:
                    # Desperate fallback: take last 2 columns? 
                    # Usually cell_id is first.
                    x = tbl.column(tbl.num_columns - 2).to_numpy(zero_copy_only=False)
                    y = tbl.column(tbl.num_columns - 1).to_numpy(zero_copy_only=False)
                else:
                    raise ValueError(f"Embedding {embedding_name} has fewer than 2 columns")

            # Cache it if small enough? Embeddings are usually small.
            # Let's cache it.
            self.embeddings_cache[embedding_name] = pd.DataFrame({col_x: x, col_y: y})
            
            return x, y
        except Exception as e:
            # Try alternative column names (sometimes just 0, 1 or without prefix?)
            # But convert_to_parquet.py is strict.
            raise ValueError(f"Embedding {embedding_name} not found at {path}: {e}")
    
    def genes(self) -> List[str]:
        """Get list of all genes"""
        # Prefer display names if we have a mapping; otherwise return internal names.
        genes = list(self.gene_to_chunk.keys())
        if self.internal_to_display:
            return [self.internal_to_display.get(g, g) for g in genes]
        return genes
    
    def obs_count(self) -> int:
        """Get number of cells"""
        return self.n_cells
    
    def var_count(self) -> int:
        """Get number of genes"""
        return self.n_genes


# Global backend instances
_backends: Dict[str, ParquetBackend] = {}


def get_backend(dataset_id: str = None) -> ParquetBackend:
    """Get or create the ParquetBackend instance for the given dataset_id"""
    global _backends
    
    # Import here to avoid circular imports
    from .config import get_datasets

    # If dataset_id is not provided, try to find a default
    if dataset_id is None:
        datasets = get_datasets()
        if datasets:
            dataset_id = datasets[0].id
        else:
            # Fallback to single instance mode using env var
            if "default" not in _backends:
                _backends["default"] = ParquetBackend()
            return _backends["default"]

    if dataset_id in _backends:
        return _backends[dataset_id]

    # Find config for this dataset
    datasets = get_datasets()
    cfg = next((d for d in datasets if d.id == dataset_id), None)
    
    if not cfg:
        # Fallback: if we can't find it in config, but it's the "default" one
        # we might want to fallback to PARQUET_DATA_DIR.
        # But let's be strict for now to ensure config is correct.
        raise ValueError(f"Dataset {dataset_id} not found in configuration")

    # Initialize backend with URI from config
    print(f"Initializing backend for {dataset_id} at {cfg.uri}")
    backend = ParquetBackend(parquet_dir=cfg.uri)
    _backends[dataset_id] = backend
    return backend
