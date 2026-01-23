# Vector Fitting Workflow for Impedance Boundary Conditions

Complete workflow for developing patient-specific `vectorFittingImpedance` boundary conditions from clinical or simulation data.

## Overview

The vector fitting workflow allows you to:
1. **Extract impedance** from CFD results or use clinical measurements
2. **Fit rational function** Z(s) = d + Σrᵢ/(s-pᵢ) to impedance data
3. **Generate OpenFOAM BC** with optimized parameters for your case

**Advantages over standard Windkessel:**
- Patient-specific tuning from clinical data
- Multi-harmonic accuracy (captures multiple resonance peaks)
- Better waveform match: 10-20% improvement for complex pulsatile flows
- Automatic parameter estimation (no manual tuning)

## Quick Start

### Method 1: From Existing CFD Simulation

```bash
# Step 1: Extract impedance from your OpenFOAM case
python utilities/extractImpedanceFromCFD.py \
    -case /home/mchi4jw4/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_implicit \
    -outlet outlet1 \
    -start-time 0.5

# This creates: impedance_outlet1.csv

# Step 2: Fit vector model (order 4 recommended)
python utilities/impedanceVectorFit.py \
    -input impedance_outlet1.csv \
    -order 4 \
    -output outlet1_BC.txt

# Step 3: Copy BC entry to your case's 0/p file
cat outlet1_BC.txt
```

### Method 2: From Clinical Impedance Data

If you have impedance measurements (e.g., from 4D Flow MRI):

```bash
# Prepare CSV with format:
# frequency[Hz], magnitude[Pa·s/m³], phase[rad]
# OR
# frequency[Hz], real[Pa·s/m³], imag[Pa·s/m³]

python utilities/impedanceVectorFit.py \
    -input clinical_impedance_outlet1.csv \
    -order 4 \
    -outlet-name outlet1 \
    -rho 1060 \
    -p0 13332
```

## Detailed Workflow

### Step 1: Impedance Data Acquisition

You have three options for impedance data:

#### Option A: Extract from CFD (Your Case)

Use your existing simulation at `/home/mchi4jw4/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_implicit`:

```bash
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure

# Extract impedance for all outlets
python utilities/extractImpedanceFromCFD.py \
    -case /home/mchi4jw4/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_implicit \
    -outlet outlet1 outlet2 outlet3 outlet4 \
    -start-time 0.5 \
    -rho 1060

# Output files:
#   impedance_outlet1.csv
#   impedance_outlet2.csv
#   impedance_outlet3.csv
#   impedance_outlet4.csv
#   timeseries_*.png (pressure/flowrate plots)
#   impedance_*.png (Bode diagrams)
```

**Parameters:**
- `-case`: Path to OpenFOAM case directory
- `-outlet`: Outlet patch name(s) from 0/p file
- `-start-time`: Exclude initial transients (recommended: 0.5-1.0s)
- `-rho`: Blood density [kg/m³], default 1060

**What it does:**
1. Reads pressure P(t) and flowrate Q(t) from all time directories
2. Applies FFT: P(t) → P̂(ω), Q(t) → Q̂(ω)
3. Computes impedance: Z(ω) = P̂(ω) / Q̂(ω)
4. Saves to CSV with frequency, magnitude, phase

**Tip:** Use `-start-time` to exclude startup transients. For cardiac simulations, exclude the first 1-2 cardiac cycles.

#### Option B: Use Clinical Impedance Measurements

If you have 4D Flow MRI or catheter-based impedance measurements:

**CSV Format 1** (magnitude-phase):
```csv
frequency,magnitude,phase
0.5,120000000,0.0
1.0,115000000,0.1
1.5,110000000,0.15
2.0,108000000,0.2
```

**CSV Format 2** (real-imaginary):
```csv
frequency,real,imag
0.5,120000000,0
1.0,114885607,11487513
1.5,106743890,28505379
2.0,101552544,41816363
```

**Units:**
- Frequency: Hz
- Magnitude/Real/Imag: Pa·s/m³ (dynamic impedance)
- Phase: radians

#### Option C: Use 1D Model Output

If you have 1D arterial network simulation output, export impedance at outlet locations.

### Step 2: Vector Fitting

Fit rational function to impedance data:

```bash
python utilities/impedanceVectorFit.py \
    -input impedance_outlet1.csv \
    -order 4 \
    -output outlet1_BC.txt \
    -outlet-name outlet1 \
    -rho 1060 \
    -p0 13332
```

**Parameters:**
- `-input`: CSV file with impedance data
- `-order`: Number of poles (recommended 4, range 4-6)
- `-output`: Output file for OpenFOAM BC entry
- `-outlet-name`: Patch name in OpenFOAM mesh
- `-rho`: Fluid density [kg/m³], default 1060
- `-p0`: Reference pressure [Pa], default 13332 (≈100 mmHg)

**Order Selection Guidelines:**

| Order | Description | When to Use |
|-------|-------------|-------------|
| 3 | Minimal | Similar to standard Windkessel, limited improvement |
| **4** | **Recommended** | Good balance: captures 1-2 harmonics, stable, fast |
| 5-6 | Advanced | Multi-harmonic cases, may require more careful tuning |
| 7+ | Specialist | Research only, risk of overfitting |

**Output Files:**
- `outlet1_BC.txt`: OpenFOAM dictionary entry (copy to 0/p)
- `outlet1_BC_bode.png`: Fitted vs. original impedance (magnitude & phase)
- `outlet1_BC_poles.png`: Pole-residue diagram

**Example Output:**
```
Vector Fitting Results
======================

Direct term d = 2.856321e+08 Pa·s/m³

Poles [rad/s]:
  p1 = -3.141593e+00
  p2 = -9.424778e+00
  p3 = -2.827433e+01
  p4 = -8.482300e+01

Residues [Pa/m³]:
  r1 = -4.523100e+07
  r2 = 3.214567e+07
  r3 = -1.892345e+07
  r4 = 8.234561e+06

Fit quality:
  RMSE: 2.34e+06 Pa·s/m³ (1.95%)
  Max error: 5.67e+06 Pa·s/m³
```

### Step 3: Validate Fitting Quality

**Check the Bode plot** (`*_bode.png`):
- Fitted curve should closely follow original data
- Pay attention to low-frequency regime (physiologically important)
- RMSE < 5% is excellent, < 10% is good

**Check the pole diagram** (`*_poles.png`):
- All poles must be in left half-plane (negative real part)
- Poles spread across frequency range
- If all poles cluster, try different order

**Common Issues:**

| Symptom | Cause | Fix |
|---------|-------|-----|
| Poor fit at low frequencies | Order too low | Increase to 5-6 |
| Oscillations in fit | Overfitting | Decrease order to 3-4 |
| Unstable poles (warning) | Noisy data | Filter input data, increase order |
| All poles identical | Bad initial guess | Add frequency points |

### Step 4: Use in OpenFOAM

Copy the BC entry from `outlet1_BC.txt` to your `0/p` file:

**Before** (standard Windkessel):
```cpp
outlet1
{
    type            modularWKPressure;
    phi             phi;
    order           3;
    R               1.2e8;      // [Pa·s/m³]
    C               1.5e-9;     // [m³/Pa]
    Z               2.4e7;      // [Pa·s/m³]
    p0              13332;      // [Pa]
    value           uniform 12.58;
}
```

**After** (vector fitting):
```cpp
outlet1
{
    type                  vectorFittingImpedance;
    phi                   phi;
    U                     U;
    couplingMode          implicit;  // Use implicit for stability
    order                 4;

    // Vector fitting parameters (DYNAMIC units - auto-converted)
    directTerm            2.8563210000e+08;  // [Pa·s/m³]
    poles                 (-3.1415930000e+00 -9.4247780000e+00 -2.8274330000e+01 -8.4823000000e+01);  // [rad/s]
    residues              (-4.5231000000e+07 3.2145670000e+07 -1.8923450000e+07 8.2345610000e+06);  // [Pa/m³]

    // Fluid properties
    rho                   1060;  // [kg/m³]

    // Initial conditions
    value                 uniform 12.577359;  // [m²/s²] kinematic = 13332/rho
}
```

**Important:**
- Change `type` from `modularWKPressure` to `vectorFittingImpedance`
- Keep `phi` name consistent
- Use `couplingMode implicit` for better stability
- Copy exact pole/residue values (high precision)

### Step 5: Run Simulation

```bash
# Recompile BC library (if not already done)
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure
wmake libso

# Run simulation
cd /path/to/your/case
pimpleFoam  # or your solver

# Monitor convergence
tail -f log.pimpleFoam
```

**Expected behavior:**
- Initial pressure may differ slightly from Windkessel
- Convergence should be similar or better (especially with implicit coupling)
- Pressure waveform will better match target impedance spectrum

## Example: Complete Workflow for BPM120 Case

Let's develop vector fitting BC for your case at `/home/mchi4jw4/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246`:

```bash
# 1. Navigate to utilities
cd /home/mchi4jw4/OpenFOAM/mchi4jw4-12/src/modularWKPressure

# 2. Extract impedance from implicit case (exclude first 0.5s transients)
python utilities/extractImpedanceFromCFD.py \
    -case /home/mchi4jw4/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_implicit \
    -outlet outlet1 outlet2 outlet3 outlet4 \
    -start-time 0.5 \
    -rho 1060

# Output:
#   impedance_outlet1.csv, impedance_outlet2.csv, impedance_outlet3.csv, impedance_outlet4.csv
#   timeseries_*.png, impedance_*.png

# 3. Fit vector models for each outlet
for outlet in outlet1 outlet2 outlet3 outlet4; do
    python utilities/impedanceVectorFit.py \
        -input impedance_${outlet}.csv \
        -order 4 \
        -output ${outlet}_vectorFitting_BC.txt \
        -outlet-name ${outlet} \
        -rho 1060 \
        -p0 10665.76
done

# Output:
#   outlet1_vectorFitting_BC.txt, ..., outlet4_vectorFitting_BC.txt
#   *_bode.png, *_poles.png

# 4. Review fit quality
ls -lh *.png
# Open Bode plots to verify good fit

# 5. Create new case with vector fitting BCs
# (Manually copy BC entries from *_vectorFitting_BC.txt to new case's 0/p)

# 6. Run new simulation
# cd /path/to/new/case
# pimpleFoam
```

## Comparison: Windkessel vs. Vector Fitting

| Aspect | 3-Element Windkessel | Vector Fitting (order 4) |
|--------|----------------------|--------------------------|
| **Parameters** | 3 (R, C, Z) | 9 (d + 4 poles + 4 residues) |
| **Tuning** | Manual or empirical formulas | Automatic from data |
| **Frequency accuracy** | Single harmonic (fundamental) | Multi-harmonic (1-2 harmonics) |
| **Waveform match** | Good for simple flows | Better for complex pulsatile |
| **Patient-specific** | Limited | Yes (from clinical data) |
| **Computational cost** | Minimal | Similar (recursive convolution) |
| **Ease of use** | Simple | Requires fitting step |
| **Recommended for** | Prototyping, simple geometries | Production, patient-specific |

## Advanced Topics

### Using Different Data Sources

**From 4D Flow MRI:**
- Export velocity field at outlet plane
- Compute flowrate Q(t) = ∫ v·n dA
- Measure pressure P(t) from PC-MRI
- FFT to get Z(ω) = P̂(ω)/Q̂(ω)
- Use CSV format: frequency, real, imag

**From 1D Models:**
- Run 1D arterial network simulation
- Export outlet impedance Z(ω) at frequency points
- Convert to CSV format
- Apply vector fitting

**From Literature Data:**
- Use published impedance spectra (e.g., Westerhof et al.)
- Digitize plots if needed
- Scale to patient-specific anatomy (use Murray's law)

### Handling Noisy Data

If impedance data is noisy:

```python
# Add to extractImpedanceFromCFD.py workflow:
from scipy.signal import savgol_filter

# Smooth pressure/flowrate before FFT
pressures = savgol_filter(pressures, window_length=11, polyorder=3)
flowrates = savgol_filter(flowrates, window_length=11, polyorder=3)
```

Or increase vector fitting order to average out noise.

### Multi-Outlet Optimization

For cases with multiple outlets, you can:
1. Fit each outlet independently (simpler)
2. Coupled fitting with flow distribution constraint (advanced)

Example for independent fitting:
```bash
# Process all outlets in parallel
for outlet in outlet1 outlet2 outlet3 outlet4; do
    (python utilities/impedanceVectorFit.py -input impedance_${outlet}.csv -order 4 &)
done
wait
```

## Troubleshooting

### "No data extracted for outlet"
- Check outlet name matches exactly in 0/p file
- Verify case has pressure `p` and flux `phi` fields
- Check file format (ASCII vs binary)

### "Unstable poles detected"
- Usually auto-corrected (forced negative)
- If persists, try different order (±1)
- Check input data for NaN or inf values

### Poor fit quality (RMSE > 10%)
- Increase order: 4 → 5 or 6
- Check input data quality (smooth time series?)
- Try different FFT window function
- Exclude more transients (increase `-start-time`)

### Simulation diverges with vector fitting BC
- Start with `couplingMode explicit` (more stable)
- Reduce timestep initially
- Switch to `couplingMode implicit` after stable
- Check pole values (very negative poles → stiff ODE)

### BC doesn't match expected waveform
- Verify impedance extraction used correct time range
- Check units in input CSV (Pa·s/m³ for impedance)
- Compare Bode plot: does fit capture key features?
- Validate reference pressure p0 is correct

## References

### Vector Fitting Algorithm
- Gustavsen, B., & Semlyen, A. (1999). "Rational approximation of frequency domain responses by vector fitting." *IEEE Transactions on Power Delivery*, 14(3), 1052-1061.

### Cardiovascular Application
- Fevola, E., Ballarin, F., Jiménez-Juan, L., Fremes, S. E., & Grivet-Talocia, S. (2023). "A vector fitting approach for the automated estimation of lumped boundary conditions of 1D circulation models." *Frontiers in Physiology*, 14, 1250204.
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC10465662/

### Impedance Measurements
- Westerhof, N., Lankhaar, J. W., & Westerhof, B. E. (2009). "The arterial Windkessel." *Medical & Biological Engineering & Computing*, 47(2), 131-141.

### 4D Flow MRI
- Markl, M., et al. (2012). "4D flow MRI." *Journal of Magnetic Resonance Imaging*, 36(5), 1015-1036.

## Support

For questions or issues:
1. Check this README and main `doc/BC_SUMMARY.md`
2. Review example plots (*.png) for fit quality
3. Validate input data format (CSV headers)
4. Check OpenFOAM compilation logs (wmake output)

Common mistakes:
- Wrong units in CSV (must be Pa·s/m³ for impedance)
- Outlet name mismatch between -outlet flag and 0/p file
- Not excluding transients (-start-time too small)
- Using binary OpenFOAM fields (extraction script needs ASCII)

---

**Version:** 1.0
**Last Updated:** 2025-11-15
**Maintainer:** OpenFOAM-WK Development Team
