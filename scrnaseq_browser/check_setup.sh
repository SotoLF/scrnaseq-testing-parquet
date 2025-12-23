#!/bin/bash
# Quick setup verification script

echo "🔍 Verificando configuración de scRNA-seq Browser..."
echo ""

# Check conda environment
echo "1️⃣ Ambiente Conda:"
if conda env list | grep -q "scrnaseq_browser"; then
    echo "   ✅ Ambiente 'scrnaseq_browser' encontrado"
else
    echo "   ❌ Ambiente 'scrnaseq_browser' NO encontrado"
    echo "   → Ejecutar: conda env create -f environment_minimal.yml"
fi
echo ""

# Check .env file
echo "2️⃣ Archivo .env:"
if [ -f .env ]; then
    echo "   ✅ Archivo .env existe"
    if grep -q "PARQUET_DATA_DIR=./parquet_data" .env; then
        echo "   ✅ Configurado para datos locales"
    elif grep -q "PARQUET_DATA_DIR=s3://" .env; then
        echo "   ℹ️  Configurado para S3"
        if grep -q "AWS_ACCESS_KEY_ID=YOUR_AWS" .env; then
            echo "   ⚠️  AWS credentials no configuradas"
        fi
    fi
else
    echo "   ❌ Archivo .env NO existe"
    echo "   → Ejecutar: cp .env.example .env"
fi
echo ""

# Check backend files
echo "3️⃣ Backend:"
if [ -f backend/app/main.py ] && [ -f backend/app/parquet_backend.py ]; then
    echo "   ✅ Archivos backend OK"
else
    echo "   ❌ Archivos backend faltantes"
fi
echo ""

# Check shiny files
echo "4️⃣ Shiny App:"
if [ -f shiny/app.R ]; then
    echo "   ✅ Archivos Shiny OK"
else
    echo "   ❌ Archivos Shiny faltantes"
fi
echo ""

# Check parquet data
echo "5️⃣ Datos Parquet:"
if [ -d parquet_data/expression ] && [ -f parquet_data/manifest.json ]; then
    chunk_count=$(ls parquet_data/expression/*.parquet 2>/dev/null | wc -l)
    echo "   ✅ Datos locales: $chunk_count chunks encontrados"
else
    echo "   ⚠️  Sin datos locales (debe usar S3)"
fi
echo ""

# Check Python packages
echo "6️⃣ Dependencias Python:"
if conda list -n scrnaseq_browser | grep -q "fastapi"; then
    echo "   ✅ FastAPI instalado"
else
    echo "   ❌ Dependencias faltantes"
    echo "   → Reinstalar: conda env create -f environment_minimal.yml --force"
fi
echo ""

# Summary
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "📋 PARA INICIAR:"
echo ""
echo "   1. conda activate scrnaseq_browser"
echo "   2. ./run_app.sh"
echo ""
echo "   Shiny: http://localhost:3838"
echo "   API:   http://localhost:8000"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

