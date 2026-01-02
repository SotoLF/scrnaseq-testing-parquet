# scRNA-seq Browser

Aplicación interactiva para visualizar y analizar datos de scRNA-seq con backend FastAPI + Parquet y frontend Shiny.

## 🚀 Características

- **Visualización UMAP/PCA** interactiva con miles de células
- **Violin plots** por gen y grupo celular
- **Backend ultra-rápido** con Parquet (~10-50ms)
- **Soporte multi-dataset**: parquet_out y parquet_bm
- **Datos locales o S3** - configurable
- **Variables categóricas**: Annotation, Cluster, Sample, etc.

## ⚡ Inicio Rápido

```bash
# 1. Activar ambiente conda
conda activate scrnaseq_browser

# 2. Configurar dataset en .env
#    PARQUET_DATA_DIR=s3://luis-scrnaseq-data/parquet_bm

# 3. Ejecutar aplicación completa
./run_app.sh

# 4. Abrir en navegador
# Shiny: http://localhost:3838
# API:   http://localhost:8000
```

## 📦 Instalación (Primera vez)

```bash
# 1. Crear ambiente conda
conda env create -f environment.yml

# 2. Activar ambiente
conda activate scrnaseq_browser

# 3. Configurar .env con tu dataset
cp .env.example .env
# Editar .env con tu configuración

# 4. Listo para usar
./run_app.sh
```

## 📊 Datasets Soportados

### parquet_out (strict_epdsc)
- **Células**: ~8,000
- **Genes**: ~8,000
- **Layers**: raw, cellbender, MAGIC_imputed_data
- **Formato**: parquet_wide_gene_chunks

### parquet_bm (BoneMarrow)
- **Células**: 75,386
- **Genes**: 26,708
- **Layers**: X, raw_counts, unlogged_normalized, lognorm_pseudocount.1
- **Variables categóricas**: Annotation, Cluster_Coarse, Cluster_Fine, Sample
- **Formato**: sparse_long (X), parquet_wide (otros)

## 🎯 Uso de la Aplicación

1. **Seleccionar gene**: Buscar en el dropdown
2. **Seleccionar color**: Variable categórica o expresión
3. **Visualizar UMAP**: Células coloreadas por grupo o expresión
4. **Violin plot**: Distribución por grupos celulares
5. **Interactivo**: Click en UMAP para seleccionar células

## �� Estructura del Proyecto

```
scrnaseq_browser/
├── backend/                    # FastAPI backend
│   ├── app/
│   │   ├── main.py            # API endpoints
│   │   ├── parquet_backend.py # Lector Parquet multi-formato
│   │   ├── cache.py           # Sistema de caché
│   │   └── config.py          # Configuración
│   └── requirements.txt
├── shiny/                     # Shiny R frontend
│   ├── app.R                  # UI y servidor
│   ├── run_app.R
│   └── www/                   # JS personalizado
├── parquet_bm/                # Dataset BoneMarrow (local)
├── parquet_out/               # Dataset strict_epdsc (local)
├── logs/                      # Logs de ejecución
├── .env                       # Configuración
├── environment.yml            # Conda environment
├── run_app.sh                 # 🚀 Script principal
└── README.md
```

## ⚙️ Configuración (.env)

```bash
# Puertos
API_PORT=8000
SHINY_PORT=3838

# Dataset - elegir uno:
# Local:
PARQUET_DATA_DIR=./parquet_bm

# O desde S3:
PARQUET_DATA_DIR=s3://luis-scrnaseq-data/parquet_bm

# AWS (si usas S3)
AWS_ACCESS_KEY_ID=tu_key
AWS_SECRET_ACCESS_KEY=tu_secret
AWS_DEFAULT_REGION=us-east-1

# Cache (opcional)
USE_REDIS=false
```

## 🔧 Scripts de Conversión

### Convertir h5ad a Parquet (dataset pequeño)
```bash
python convert_to_parquet.py --input data.h5ad --output parquet_out
```

### Convertir h5ad grande con streaming (bajo RAM)
```bash
python h5ad_to_parquet_stream_h5py.py \
  --h5ad BoneMarrow.h5ad \
  --outdir parquet_bm \
  --matrices X layers/raw_counts layers/unlogged_normalized \
  --obsm X_umap X_pca \
  --chunk-genes 256
```

### Agregar variables categóricas a parquet existente
```bash
python add_categorical_to_obs.py \
  --h5ad BoneMarrow.h5ad \
  --parquet-dir parquet_bm \
  --backup
```

## 🛠️ Troubleshooting

**Error: conda environment not found**
```bash
conda env create -f environment.yml
conda activate scrnaseq_browser
```

**Backend no responde**
```bash
# Ver logs
tail -f logs/backend.log

# Verificar configuración
./check_setup.sh
```

**Puerto ya en uso**
```bash
# Cambiar puerto en .env
API_PORT=8001
SHINY_PORT=3839

# O limpiar puertos
CLEAN_START=1 ./run_app.sh
```

**Ver logs en tiempo real**
```bash
tail -f logs/backend.log
tail -f logs/shiny.log
```

## 📊 Performance

| Operación | Tiempo |
|-----------|--------|
| Gene retrieval | ~10-50ms |
| UMAP coordinates | ~20ms |
| Cell metadata | ~30ms |

- **256 genes/chunk** para balance memoria/velocidad
- **Cache en memoria** para queries repetidas
- **Streaming desde S3** para datasets grandes

## 🔗 API Endpoints

| Endpoint | Descripción |
|----------|-------------|
| `GET /genes` | Lista de genes disponibles |
| `GET /obs_columns` | Columnas de metadatos |
| `GET /obsm_keys` | Embeddings disponibles (UMAP, PCA) |
| `GET /obsm_xy?key=X_umap` | Coordenadas de embedding |
| `GET /obs_col?col=Annotation` | Valores de columna |
| `POST /expr_batch` | Expresión de múltiples genes |
| `GET /layers` | Layers disponibles |

Documentación interactiva: http://localhost:8000/docs

## 📝 Licencia

MIT License
