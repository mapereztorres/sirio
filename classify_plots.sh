#!/bin/bash

# Crear carpetas destino y subcarpetas por geometría - PDF
mkdir -p OUTPUT/T_tauri_proposal
mkdir -p OUTPUT/T_tauri_proposal/closed-dipole
mkdir -p OUTPUT/T_tauri_proposal/open-parker-spiral
mkdir -p OUTPUT/T_tauri_proposal/pfss

# Crear carpetas destino y subcarpetas por geometría - PNG
mkdir -p OUTPUT/T_tauri_proposal_png
mkdir -p OUTPUT/T_tauri_proposal_png/closed-dipole-png
mkdir -p OUTPUT/T_tauri_proposal_png/open-parker-spiral-png
mkdir -p OUTPUT/T_tauri_proposal_png/pfss-png

# Lista de directorios con los nombres reales
dirs=(
    "2MASS_J02281782+7238009_Earth"
    "2MASS_J02281782+7238009_Jupiter"
    "2MASS_J04244506+2701447_Earth"
    "2MASS_J04244506+2701447_Jupiter"
    "2MASS_J22361978+4006273_Earth"
    "2MASS_J22361978+4006273_Jupiter"
    "HD_283447_Earth"
    "HD_283447_Jupiter"
    "IC_348_LRL_224_Earth"
    "IC_348_LRL_224_Jupiter"
    "LkCa_4_Earth"
    "LkCa_4_Jupiter"
    "T_Tau_Earth"
    "T_Tau_Jupiter"
    "UCAC4_591-011866_Earth"
    "UCAC4_591-011866_Jupiter"
    "V1320_Tau_Earth"
    "V1320_Tau_Jupiter"
)

# Función para clasificar archivos por geometría
classify() {
    local file="$1"
    local dest_base="$2"
    local suffix="$3"
    local filename
    filename=$(basename "$file")

    # Copiar a la carpeta general
    cp "$file" "$dest_base/$filename"

    # Copiar por geometría
    if [[ "$filename" == *"closed-dipole"* ]]; then
        cp "$file" "$dest_base/closed-dipole${suffix}/$filename"
    elif [[ "$filename" == *"open-parker-spiral"* ]]; then
        cp "$file" "$dest_base/open-parker-spiral${suffix}/$filename"
    elif [[ "$filename" == *"pfss"* ]]; then
        cp "$file" "$dest_base/pfss${suffix}/$filename"
    else
        echo "⚠️  Geometría no reconocida: $filename"
    fi
}

# Procesar cada directorio
for dir in "${dirs[@]}"; do
    full_dir="OUTPUT/$dir"

    if [ ! -d "$full_dir" ]; then
        echo "⚠️  No encontrado: $full_dir"
        continue
    fi

    # --- PDFs ---
    if [ -d "$full_dir/FLUX_PDF" ]; then
        for pdf in "$full_dir/FLUX_PDF"/*.pdf; do
            [ -f "$pdf" ] || continue
            classify "$pdf" "OUTPUT/T_tauri_proposal" ""
        done
        echo "✅ $dir — PDFs procesados"
    else
        echo "⚠️  Sin carpeta FLUX_PDF en: $full_dir"
    fi

    # --- PNGs ---
    if [ -d "$full_dir/FLUX_PNG" ]; then
        for png in "$full_dir/FLUX_PNG"/*.png; do
            [ -f "$png" ] || continue
            classify "$png" "OUTPUT/T_tauri_proposal_png" "-png"
        done
        echo "✅ $dir — PNGs procesados"
    else
        echo "⚠️  Sin carpeta FLUX_PNG en: $full_dir"
    fi

done

# Copiar latex_table.tex
if [ -f "OUTPUT/latex_table.tex" ]; then
    cp OUTPUT/latex_table.tex OUTPUT/T_tauri_proposal/
    echo "✅ latex_table.tex copiado"
else
    echo "⚠️  No encontrado: OUTPUT/latex_table.tex"
fi

echo ""
echo "✅ Listo. Estructura final:"
echo "   OUTPUT/T_tauri_proposal/                     ← todos los PDFs + latex_table.tex"
echo "   OUTPUT/T_tauri_proposal/closed-dipole/"
echo "   OUTPUT/T_tauri_proposal/open-parker-spiral/"
echo "   OUTPUT/T_tauri_proposal/pfss/"
echo ""
echo "   OUTPUT/T_tauri_proposal_png/                 ← todos los PNGs"
echo "   OUTPUT/T_tauri_proposal_png/closed-dipole-png/"
echo "   OUTPUT/T_tauri_proposal_png/open-parker-spiral-png/"
echo "   OUTPUT/T_tauri_proposal_png/pfss-png/"
