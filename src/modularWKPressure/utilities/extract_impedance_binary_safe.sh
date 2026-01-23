#!/bin/bash
#
# Extract Impedance from Binary OpenFOAM Case
# ============================================
#
# This script works with binary format OpenFOAM cases by using
# OpenFOAM's native postProcess tools, which can read binary files.
#
# Usage:
#   bash extract_impedance_binary_safe.sh <case_directory> [start_time]
#
# Example:
#   bash extract_impedance_binary_safe.sh ~/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_stable 0.5

set -e

CASE_DIR="$1"
START_TIME="${2:-0.5}"
OUTLETS="${3:-outlet1 outlet2 outlet3 outlet4}"
RHO="${4:-1060}"

if [ -z "$CASE_DIR" ]; then
    echo "Usage: $0 <case_directory> [start_time] [outlets] [rho]"
    echo ""
    echo "Example:"
    echo "  $0 ~/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_stable 0.5"
    echo ""
    exit 1
fi

if [ ! -d "$CASE_DIR" ]; then
    echo "Error: Case directory not found: $CASE_DIR"
    exit 1
fi

echo "========================================================================"
echo "Extract Impedance from Binary OpenFOAM Case"
echo "========================================================================"
echo ""
echo "Case:       $CASE_DIR"
echo "Start time: $START_TIME s"
echo "Outlets:    $OUTLETS"
echo "Density:    $RHO kg/m³"
echo ""

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ============================================================================
# STEP 1: Create functionObjects for extraction
# ============================================================================

echo "========================================================================"
echo "STEP 1: Creating functionObjects for extraction"
echo "========================================================================"
echo ""

cd "$CASE_DIR"

# Create extraction function dictionary (OpenFOAM 12 format)
cat > system/extractOutletData << 'EOF'
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Version:  12
     \\/     M anipulation  |
\*---------------------------------------------------------------------------*/
FoamFile
{
    format      ascii;
    class       dictionary;
    location    "system";
    object      extractOutletData;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

outlet1Pressure
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (p);
    operation       average;
    regionType      patch;
    name            outlet1;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet1Flow
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (phi);
    operation       sum;
    regionType      patch;
    name            outlet1;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet2Pressure
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (p);
    operation       average;
    regionType      patch;
    name            outlet2;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet2Flow
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (phi);
    operation       sum;
    regionType      patch;
    name            outlet2;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet3Pressure
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (p);
    operation       average;
    regionType      patch;
    name            outlet3;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet3Flow
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (phi);
    operation       sum;
    regionType      patch;
    name            outlet3;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet4Pressure
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (p);
    operation       average;
    regionType      patch;
    name            outlet4;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

outlet4Flow
{
    type            surfaceFieldValue;
    libs            ("libfieldFunctionObjects.so");

    fields          (phi);
    operation       sum;
    regionType      patch;
    name            outlet4;

    writeControl    timeStep;
    writeInterval   1;
    writeFields     false;
}

// ************************************************************************* //
EOF

echo "✓ Created system/extractOutletData"
echo ""

# ============================================================================
# STEP 2: Run postProcess to extract data
# ============================================================================

echo "========================================================================"
echo "STEP 2: Running postProcess (handles binary format)"
echo "========================================================================"
echo ""

echo "Extracting data (this may take a minute)..."

# OpenFOAM 12: Use foamPostProcess with -dict option
if command -v foamPostProcess &> /dev/null; then
    foamPostProcess -dict system/extractOutletData > postProcess.log 2>&1
elif command -v postProcess &> /dev/null; then
    postProcess -dict system/extractOutletData > postProcess.log 2>&1
else
    echo "ERROR: Neither foamPostProcess nor postProcess found!"
    exit 1
fi

if [ $? -ne 0 ]; then
    echo "ERROR: postProcess failed!"
    echo ""
    cat postProcess.log
    exit 1
fi

echo "✓ Data extracted to postProcessing/"
echo ""

# Check output
if [ -d "postProcessing" ]; then
    echo "Generated directories:"
    ls -d postProcessing/outlet* 2>/dev/null || echo "  (none found)"
    echo ""
fi

cd "$SCRIPT_DIR"

# ============================================================================
# STEP 3: Convert postProcessing data to impedance CSV
# ============================================================================

echo "========================================================================"
echo "STEP 3: Converting to impedance CSV format"
echo "========================================================================"
echo ""

# Check if Python script exists
if [ ! -f "convert_postProcessing_to_impedance.py" ]; then
    echo "ERROR: convert_postProcessing_to_impedance.py not found!"
    echo "       Expected at: $SCRIPT_DIR/convert_postProcessing_to_impedance.py"
    exit 1
fi

# Activate Python environment if available
if [ -d "/home/mchi4jw4/GitHub/AortaCFD-app/venv" ]; then
    source /home/mchi4jw4/GitHub/AortaCFD-app/venv/bin/activate
fi

# Run conversion
python3 convert_postProcessing_to_impedance.py \
    -case "$CASE_DIR" \
    -outlet $OUTLETS \
    -start-time $START_TIME \
    -rho $RHO

if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: Conversion failed!"
    exit 1
fi

echo ""
echo "========================================================================"
echo "EXTRACTION COMPLETE!"
echo "========================================================================"
echo ""
echo "Generated files:"
for outlet in $OUTLETS; do
    if [ -f "impedance_${outlet}.csv" ]; then
        echo "  ✓ impedance_${outlet}.csv"
    else
        echo "  ✗ impedance_${outlet}.csv (missing!)"
    fi
done
echo ""
echo "Next step: Run vector fitting"
echo "  python impedanceVectorFit.py -input impedance_outlet1.csv -order 4 ..."
echo ""
echo "Or use the automated workflow:"
echo "  bash run_vector_fitting_workflow.sh"
echo "  (make sure to update it to skip extraction step)"
echo ""
