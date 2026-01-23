#!/bin/bash
#
# Extract Outlet Data Using OpenFOAM postProcess
# ===============================================
#
# This uses OpenFOAM's native tools to extract data (handles binary format)

set -e

CASE_DIR="$1"
START_TIME="${2:-0.5}"

if [ -z "$CASE_DIR" ]; then
    echo "Usage: $0 <case_directory> [start_time]"
    exit 1
fi

cd "$CASE_DIR"

echo "========================================================================"
echo "Extracting outlet data using OpenFOAM postProcess"
echo "========================================================================"
echo ""
echo "Case: $CASE_DIR"
echo "Start time: $START_TIME"
echo ""

# Create functionObjects to extract patch data
echo "Creating extraction functions..."

cat > system/extractOutletData << 'EOF'
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
EOF

# Run postProcess for all time steps
echo "Running postProcess (this may take a minute)..."
postProcess -func extractOutletData > postProcess.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: postProcess failed!"
    cat postProcess.log
    exit 1
fi

echo "✓ Data extracted to postProcessing/"
echo ""

# Check output
if [ -d "postProcessing" ]; then
    echo "Generated files:"
    ls -lh postProcessing/
    echo ""
    echo "Data locations:"
    find postProcessing -name "*.dat" | head -10
fi

echo ""
echo "Next step: Convert postProcessing data to CSV for Python script"
