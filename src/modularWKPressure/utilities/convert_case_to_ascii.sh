#!/bin/bash
#
# Convert OpenFOAM Case from Binary to ASCII Format
# ==================================================
#
# The Python impedance extraction script requires ASCII format files.
# This script converts a case to ASCII format for processing.
#
# Usage:
#   bash convert_case_to_ascii.sh <case_directory>

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <case_directory>"
    echo ""
    echo "Example:"
    echo "  $0 /path/to/openfoam_stable"
    exit 1
fi

CASE_DIR="$1"

if [ ! -d "$CASE_DIR" ]; then
    echo "Error: Directory not found: $CASE_DIR"
    exit 1
fi

echo "========================================================================"
echo "Converting OpenFOAM Case to ASCII Format"
echo "========================================================================"
echo ""
echo "Case: $CASE_DIR"
echo ""

# Check if case is already ASCII
FORMAT=$(grep "format" "$CASE_DIR/0/p" 2>/dev/null | head -1 | awk '{print $2}' | tr -d ';')
echo "Initial 0/p format: $FORMAT"

# Check time directories
TIME_DIR_FORMAT=$(find "$CASE_DIR" -maxdepth 1 -type d -name "[0-9]*" -o -name "0.[0-9]*" | head -1)
if [ ! -z "$TIME_DIR_FORMAT" ]; then
    TIME_DIR=$(basename "$TIME_DIR_FORMAT")
    if [ -f "$CASE_DIR/$TIME_DIR/p" ]; then
        TIME_FORMAT=$(grep "format" "$CASE_DIR/$TIME_DIR/p" 2>/dev/null | head -1 | awk '{print $2}' | tr -d ';')
        echo "Time directory format ($TIME_DIR): $TIME_FORMAT"
    fi
fi

if [ "$TIME_FORMAT" == "ascii" ]; then
    echo ""
    echo "✓ Case is already in ASCII format. No conversion needed."
    exit 0
fi

echo ""
echo "Converting to ASCII format..."
echo ""

# Use foamFormatConvert to convert all fields to ASCII
cd "$CASE_DIR"

if command -v foamFormatConvert &> /dev/null; then
    echo "Using foamFormatConvert..."
    foamFormatConvert -ascii
    CONVERSION_STATUS=$?
else
    echo "foamFormatConvert not found. Trying alternative method..."

    # Alternative: Use foamDictionary to convert each file
    for time_dir in $(ls -d [0-9]* 0.[0-9]* 2>/dev/null); do
        echo "  Converting $time_dir..."
        for field in p U phi; do
            if [ -f "$time_dir/$field" ]; then
                # Check if binary
                field_format=$(grep "format" "$time_dir/$field" | head -1 | awk '{print $2}' | tr -d ';')
                if [ "$field_format" == "binary" ]; then
                    # Use foamDictionary to read and write (converts to ASCII)
                    foamDictionary -entry format -set ascii "$time_dir/$field" > /dev/null 2>&1 || echo "    Warning: Could not convert $time_dir/$field"
                fi
            fi
        done
    done
    CONVERSION_STATUS=0
fi

cd - > /dev/null

if [ $CONVERSION_STATUS -eq 0 ]; then
    echo ""
    echo "✓ Conversion completed successfully"
    echo ""

    # Verify conversion
    if [ ! -z "$TIME_DIR" ]; then
        NEW_FORMAT=$(grep "format" "$CASE_DIR/$TIME_DIR/p" 2>/dev/null | head -1 | awk '{print $2}' | tr -d ';')
        echo "Verification: $TIME_DIR/p format is now: $NEW_FORMAT"
    fi

    echo ""
    echo "Note: Files are now larger (ASCII takes more space than binary)"
    echo "      You can delete time directories after impedance extraction"
else
    echo ""
    echo "Error: Conversion failed"
    exit 1
fi
