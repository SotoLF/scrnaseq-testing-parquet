# scRNA-seq Browser - Guía Rápida

## ✅ Estado del Proyecto

Proyecto completamente funcional con soporte para múltiples datasets.

## 🚀 Inicio Super Rápido

```bash
# 1. Activar ambiente conda
conda activate scrnaseq_browser

# 2. Ejecutar todo
./run_app.sh

# 3. Abrir navegador
# http://localhost:3838
```

**¡Eso es todo!** 🎉

## 📦 Primera Instalación

```bash
# Crear ambiente conda (solo una vez)
conda env create -f environment.yml

# Activar y listo
conda activate scrnaseq_browser
./run_app.sh
```

## 📊 Datasets Disponibles

### Dataset 1: parquet_out (strict_epdsc)
```bash
# En .env:
PARQUET_DATA_DIR=s3://luis-scrnaseq-data/parquet_out
# o local:
PARQUET_DATA_DIR=./parquet_out
```
- ~8,000 células
- Layers: raw, cellbender, MAGIC_imputed_data

### Dataset 2: parquet_bm (BoneMarrow) ⭐ Recomendado
```bash
# En .env:
PARQUET_DATA_DIR=s3://luis-scrnaseq-data/parquet_bm
# o local:
PARQUET_DATA_DIR=./parquet_bm
```
- 75,386 células, 26,708 genes
- Layers: X, raw_counts, unlogged_normalized, lognorm_pseudocount.1
- Variables categóricas: **Annotation**, Cluster_Coarse, Cluster_Fine, Sample

## 🎯 Funcionalidades

| Feature | Descripción |
|---------|-------------|
| 🔍 Gene Search | Búsqueda ultra-rápida (~10-50ms) |
| 📊 UMAP Plot | Visualización interactiva de células |
| 🎻 Violin Plot | Distribución de expresión por grupos |
| 🏷️ Color by Category | Colorear por Annotation, Cluster, etc. |
| 🖱️ Selection | Click en UMAP para filtrar células |

## ⚙️ Configuración Rápida

### Archivo .env
```bash
# Puertos
API_PORT=8000
SHINY_PORT=3838

# Dataset (elegir uno)
PARQUET_DATA_DIR=s3://luis-scrnaseq-data/parquet_bm

# AWS credentials (si usas S3)
AWS_ACCESS_KEY_ID=tu_key
AWS_SECRET_ACCESS_KEY=tu_secret
AWS_DEFAULT_REGION=us-east-1
```

## 📁 Estructura

```
scrnaseq_browser/
├── backend/app/           # FastAPI backend
│   ├── main.py           # Endpoints
│   └── parquet_backend.py # Lector multi-formato
├── shiny/                 # Shiny R frontend
│   └── app.R
├── parquet_bm/           # Dataset BoneMarrow
├── parquet_out/          # Dataset strict_epdsc
├── .env                  # Configuración
├── environment.yml       # Conda environment
└── run_app.sh           # 🚀 Script principal
```

## 🔧 Scripts Útiles

### Convertir nuevo dataset h5ad
```bash
# Para datasets grandes (bajo RAM):
python h5ad_to_parquet_stream_h5py.py \
  --h5ad MiDataset.h5ad \
  --outdir parquet_nuevo \
  --matrices X layers/raw_counts \
  --obsm X_umap X_pca
```

### Agregar variables categóricas faltantes
```bash
python add_categorical_to_obs.py \
  --h5ad MiDataset.h5ad \
  --parquet-dir parquet_nuevo
```

### Subir a S3
```bash
aws s3 sync parquet_nuevo s3://luis-scrnaseq-data/parquet_nuevo
```

## 🛠️ Troubleshooting

| Problema | Solución |
|----------|----------|
| Backend no responde | `tail -f logs/backend.log` |
| Puerto en uso | `CLEAN_START=1 ./run_app.sh` |
| No muestra categorías | Ejecutar `add_categorical_to_obs.py` |
| UMAP muestra cruz | Verificar obsm en manifest |

### Comandos de diagnóstico
```bash
# Verificar setup completo
./check_setup.sh

# Test solo backend
./test_backend.sh

# Ver logs
tail -f logs/backend.log
tail -f logs/shiny.log
```

## 🌐 URLs

| Servicio | URL |
|----------|-----|
| Shiny App | http://localhost:3838 |
| API | http://localhost:8000 |
| API Docs | http://localhost:8000/docs |

## 📈 Performance

- **Gene retrieval**: ~10-50ms
- **UMAP load**: ~20ms  
- **75k células**: sin problema
- **Cache**: queries repetidas instantáneas
