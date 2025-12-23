#!/bin/bash
# test_backend.sh - Quick test that backend can load parquet data

cd "$(dirname "$0")"

echo "🧪 Testing backend initialization..."
echo ""

cd backend
python3 -c "
import sys
sys.path.insert(0, '.')
from app.parquet_backend import ParquetBackend
import os

os.environ['PARQUET_DATA_DIR'] = './parquet_data'
try:
    backend = ParquetBackend()
    print('✅ Backend test PASSED')
    print(f'   Cells: {backend.n_cells:,}')
    print(f'   Genes: {backend.n_genes:,}')
    print(f'   Path: {backend.parquet_dir}')
except Exception as e:
    print(f'❌ Backend test FAILED: {e}')
    sys.exit(1)
" 2>&1 | grep -v "UserWarning"

echo ""
echo "Ready to run ./run_app.sh"
