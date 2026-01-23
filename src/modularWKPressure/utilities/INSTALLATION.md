# Installation Guide for Vector Fitting Utilities

## Python Dependencies

The vector fitting utilities require Python 3.8+ with scientific computing libraries.

### Quick Install

```bash
# Navigate to utilities directory
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure/utilities

# Install dependencies
pip3 install -r requirements.txt

# Or install individually:
pip3 install numpy scipy matplotlib pandas
```

### Using Virtual Environment (Recommended)

```bash
# Create virtual environment
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure/utilities
python3 -m venv venv

# Activate environment
source venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Test installation
python impedanceVectorFit.py --help
python extractImpedanceFromCFD.py --help
```

### Using Existing Environment

If you have an existing Python environment (e.g., from AortaCFD-app):

```bash
# Activate existing environment
source /home/mchi4jw4/GitHub/AortaCFD-app/venv/bin/activate

# Install missing dependencies if needed
pip install scipy matplotlib pandas

# Run utilities
python utilities/impedanceVectorFit.py --help
```

## Verification

Test that everything works:

```bash
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure/utilities

# Test vector fitting with example data
python impedanceVectorFit.py \
    -input example_impedance_aortic.csv \
    -order 4 \
    -output test_BC.txt

# Check output
cat test_BC.txt
ls -lh test_BC_*.png
```

Expected output:
- `test_BC.txt`: OpenFOAM BC entry
- `test_BC_bode.png`: Bode plot showing fit quality
- `test_BC_poles.png`: Pole-residue diagram

## Troubleshooting

### "ModuleNotFoundError: No module named 'pandas'"

Install pandas:
```bash
pip3 install pandas
```

### "Permission denied" when installing packages

Use `--user` flag:
```bash
pip3 install --user -r requirements.txt
```

Or use virtual environment (see above).

### "python3: command not found"

Your system might use `python` instead of `python3`:
```bash
python -m pip install -r requirements.txt
python impedanceVectorFit.py --help
```

### Dependencies conflict with system packages

Use virtual environment to isolate dependencies:
```bash
python3 -m venv venv_vectorfit
source venv_vectorfit/bin/activate
pip install -r requirements.txt
```

## System Requirements

- Python 3.8 or newer
- ~200 MB disk space for dependencies
- NumPy, SciPy, Matplotlib, Pandas (installed via requirements.txt)

## Usage After Installation

Once dependencies are installed, refer to:
- `README_VECTOR_FITTING.md` for complete workflow
- `impedanceVectorFit.py --help` for command-line options
- `extractImpedanceFromCFD.py --help` for extraction options
