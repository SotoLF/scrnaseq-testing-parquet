#!/bin/bash
# run_app.sh - Ejecutar scRNA-seq Browser completo
# Inicia backend FastAPI y Shiny app simultáneamente

set -e  # Exit on error

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_ROOT"

echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BLUE}  🧬 scRNA-seq Browser Launcher${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Check conda environment
if [[ "$CONDA_DEFAULT_ENV" != "scrnaseq_browser" ]]; then
    echo -e "${YELLOW}⚠️  Conda environment 'scrnaseq_browser' no está activado${NC}"
    echo -e "${YELLOW}   Activando ambiente...${NC}"
    
    # Try to activate conda env
    if command -v conda &> /dev/null; then
        eval "$(conda shell.bash hook)"
        conda activate scrnaseq_browser 2>/dev/null || {
            echo -e "${RED}❌ Error: No se pudo activar 'scrnaseq_browser'${NC}"
            echo -e "${RED}   Ejecuta: conda env create -f environment_minimal.yml${NC}"
            exit 1
        }
    else
        echo -e "${RED}❌ Error: conda no encontrado${NC}"
        exit 1
    fi
fi

echo -e "${GREEN}✅ Ambiente conda: $CONDA_DEFAULT_ENV${NC}"
echo ""

# Check .env file
if [ ! -f .env ]; then
    echo -e "${YELLOW}⚠️  Archivo .env no encontrado${NC}"
    echo -e "${YELLOW}   Creando desde .env.example...${NC}"
    cp .env.example .env
fi

# Source environment variables
set -a
source .env
set +a

echo -e "${GREEN}✅ Configuración cargada${NC}"
echo -e "   API Port: ${YELLOW}${API_PORT}${NC}"
echo -e "   Shiny Port: ${YELLOW}${SHINY_PORT}${NC}"
echo -e "   Parquet Data: ${YELLOW}${PARQUET_DATA_DIR}${NC}"
echo ""

# Ensure R dependencies are available in this conda env
echo -e "${BLUE}📦 Verificando dependencias R...${NC}"
set +e
MISSING_R_PKGS=$(Rscript -e 'pkgs <- c("shiny","bslib","httr","jsonlite","base64enc","viridisLite","ggplot2","scales","htmlwidgets","dplyr"); m <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; if(length(m)){cat(paste(m, collapse=" ")); quit(status=1)}' 2>/dev/null)
R_CHECK_STATUS=$?
set -e

if [ $R_CHECK_STATUS -ne 0 ]; then
        echo -e "${YELLOW}⚠️  Faltan paquetes R: ${MISSING_R_PKGS}${NC}"
        echo -e "${YELLOW}   Instalando via conda-forge (binarios)...${NC}"
        conda install -y -c conda-forge \
            r-shiny r-bslib r-httr r-jsonlite r-base64enc r-viridislite r-ggplot2 r-scales r-htmlwidgets r-dplyr
        echo -e "${GREEN}   ✓ Paquetes R instalados${NC}"
else
        echo -e "${GREEN}   ✓ Paquetes R OK${NC}"
fi
echo ""

# Create log directory
mkdir -p logs

# Optional: ensure a clean start by killing anything already bound to the ports.
# Usage: CLEAN_START=1 ./run_app.sh
# Default to 1 (always clean) to avoid zombie processes during development
CLEAN_START=${CLEAN_START:-1}

kill_port() {
    local port="$1"
    if command -v fuser >/dev/null 2>&1; then
        fuser -k "${port}/tcp" >/dev/null 2>&1 || true
    elif command -v lsof >/dev/null 2>&1; then
        lsof -ti "tcp:${port}" | xargs -r kill >/dev/null 2>&1 || true
    fi
}

if [ "$CLEAN_START" = "1" ]; then
    echo -e "${YELLOW}🧹 CLEAN_START=1: liberando puertos ${API_PORT} y ${SHINY_PORT}...${NC}"
    kill_port "${API_PORT}"
    kill_port "${SHINY_PORT}"
    # Also kill common orphan patterns (best-effort)
    pkill -f "R.*--file=run_app\\.R" >/dev/null 2>&1 || true
    pkill -f "uvicorn app\\.main:app" >/dev/null 2>&1 || true
    echo ""
fi

# Prefer uvicorn from the active conda env to avoid spawning processes under system Python.
PYTHON_BIN="$(command -v python)"

# Cleanup function
cleanup() {
    echo ""
    echo -e "${YELLOW}🛑 Deteniendo servicios...${NC}"
    
    # Kill background jobs
    if [ ! -z "$BACKEND_PID" ]; then
        # Kill entire process group (backend + any children)
        kill -TERM -- -$BACKEND_PID 2>/dev/null || kill $BACKEND_PID 2>/dev/null || true
        echo -e "${GREEN}   ✓ Backend detenido${NC}"
    fi
    
    if [ ! -z "$SHINY_PID" ]; then
        # Kill entire process group (shiny + any children)
        kill -TERM -- -$SHINY_PID 2>/dev/null || kill $SHINY_PID 2>/dev/null || true
        echo -e "${GREEN}   ✓ Shiny detenido${NC}"
    fi

    # Fallback (best-effort): kill orphaned instances from previous runs.
    pkill -f "R.*--file=run_app\\.R" >/dev/null 2>&1 || true
    pkill -f "uvicorn app\\.main:app" >/dev/null 2>&1 || true
    
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${GREEN}✨ Aplicación cerrada${NC}"
    exit 0
}

# Trap Ctrl+C
trap cleanup SIGINT SIGTERM

# Start backend
echo -e "${BLUE}🚀 Iniciando Backend FastAPI...${NC}"
cd backend
# Set PROJECT_ROOT for backend to resolve relative paths
# Use absolute path for logs
LOG_FILE="$PROJECT_ROOT/logs/backend.log"
PROJECT_ROOT="$PROJECT_ROOT" setsid "$PYTHON_BIN" -m uvicorn app.main:app --host 0.0.0.0 --port ${API_PORT} --log-level info > "$LOG_FILE" 2>&1 &
BACKEND_PID=$!
cd ..

echo -e "${GREEN}   ✓ Backend iniciado (PID: $BACKEND_PID)${NC}"
echo -e "   Log: logs/backend.log"
echo ""

# Wait for backend to be ready
echo -e "${YELLOW}⏳ Esperando que backend esté listo...${NC}"
sleep 3

# Check if backend is responding and actually initialized
for i in {1..10}; do
    if curl --noproxy '*' -sf http://localhost:${API_PORT}/datasets > /dev/null 2>&1; then
        echo -e "${GREEN}   ✓ Backend listo en http://localhost:${API_PORT}${NC}"
        break
    fi
    
    if [ $i -eq 10 ]; then
        echo -e "${RED}❌ Backend no responde${NC}"
        echo -e "${RED}   Revisa logs/backend.log${NC}"
        cleanup
    fi
    
    sleep 1
done

echo ""

# Start Shiny
echo -e "${BLUE}🎨 Iniciando Shiny App...${NC}"
cd shiny
SHINY_LOG="$PROJECT_ROOT/logs/shiny.log"
API_BASE_URL=http://localhost:${API_PORT} SHINY_PORT=${SHINY_PORT} setsid Rscript run_app.R > "$SHINY_LOG" 2>&1 &
SHINY_PID=$!
cd ..

echo -e "${GREEN}   ✓ Shiny iniciado (PID: $SHINY_PID)${NC}"
echo -e "   Log: logs/shiny.log"
echo ""

# Wait for Shiny to be ready
echo -e "${YELLOW}⏳ Esperando que Shiny esté listo...${NC}"
sleep 5

# Check if Shiny is responding
for i in {1..15}; do
    if curl --noproxy '*' -sf http://localhost:${SHINY_PORT}/ > /dev/null 2>&1; then
        echo -e "${GREEN}   ✓ Shiny respondiendo en http://localhost:${SHINY_PORT}${NC}"
        break
    fi

    if [ $i -eq 15 ]; then
        echo -e "${RED}❌ Shiny no responde en http://localhost:${SHINY_PORT}${NC}"
        echo -e "${RED}   Revisa logs/shiny.log${NC}"
        echo -e "${YELLOW}--- Tail logs/shiny.log ---${NC}"
        tail -n 120 logs/shiny.log || true
        cleanup
    fi

    sleep 1
done

echo ""
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${GREEN}✨ ¡Aplicación corriendo!${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""
echo -e "  🌐 Shiny App:  ${GREEN}http://localhost:${SHINY_PORT}${NC}"
echo -e "  🔧 API:        ${GREEN}http://localhost:${API_PORT}${NC}"
echo -e "  📊 API Docs:   ${GREEN}http://localhost:${API_PORT}/docs${NC}"
echo ""
echo -e "${YELLOW}📝 Logs:${NC}"
echo -e "   Backend: tail -f logs/backend.log"
echo -e "   Shiny:   tail -f logs/shiny.log"
echo ""
echo -e "${YELLOW}Press Ctrl+C to stop all services${NC}"
echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo ""

# Wait for processes
wait
