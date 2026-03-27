# CoA Test Case: Coarctation of Aorta with Windkessel Outlets

This tutorial demonstrates the 3-element Windkessel (RCR) pressure boundary condition with two-parameter backflow stabilisation on a patient-specific coarctation of aorta (CoA) geometry.

## Prerequisites

1. Compiled `libmodularWKPressure.so` (see main README)
2. OpenFOAM 12 (Foundation)
3. A pre-generated aortic mesh with patches: `inlet`, `outlet1`-`outlet4`, `wall_aorta`

**Note:** The mesh is NOT included due to size. Generate it using AortaCFD-app:

```bash
# Using AortaCFD-app:
python run_patient.py <case_id> --steps case,mesh
cp -r output/<case_id>/run_xxx/openfoam/constant/polyMesh tutorials/CoA_test/constant/

# Or provide any aortic mesh with the correct patch names.
```

## Boundary Conditions

### Pressure (`0/p`): modularWKPressure

All parameters are in **kinematic units** (divided by rho = 1060 kg/m^3):

```cpp
outlet1     // Descending aorta (largest outlet)
{
    type            modularWKPressure;
    phi             phi;
    U               U;
    couplingMode    implicit;
    order           3;              // BDF3 for accuracy
    R               268765.33;      // Distal resistance [m^-1.s^-1]
    C               3.72e-6;        // Compliance [m.s^2]
    Z               268936.63;      // Characteristic impedance [m^-1.s^-1]
    p0              10.06;          // ~80 mmHg diastolic [m^2/s^2]
    value           uniform 10.06;
}
```

Each outlet has different R/C/Z values based on Murray's law flow distribution. Smaller outlets have higher resistance and lower compliance.

### Velocity (`0/U`): stabilizedWindkesselVelocity

Two-parameter directional stabilisation:

```cpp
outlet1
{
    type                  stabilizedWindkesselVelocity;
    phi                   phi;
    betaT                 0.2;    // Tangential damping (suppresses backflow vortices)
    betaN                 0.0;    // Normal: free for Windkessel pressure-flow coupling
    enableStabilization   true;
    value                 uniform (0 0 0);
}
```

- **betaT = 0.2**: Damps tangential velocity during backflow. Increase to 0.3 for severe cases.
- **betaN = 0.0**: Keeps normal velocity free to respond to the Windkessel pressure model.

## Critical fvSolution Requirement

The Windkessel BC requires `pFinal = 1.0` and `UFinal = 1.0`:

```cpp
relaxationFactors
{
    fields    { p 0.3; pFinal 1.0; }
    equations { U 0.7; UFinal 1.0; }
}
```

Under-relaxing the final PIMPLE iteration corrupts the Windkessel pressure-flow coupling. See the main README for details.

## Running

```bash
source /opt/openfoam12/etc/bashrc
./Allrun
```

## Expected Behaviour

| Metric | Without Stabilisation | With Stabilisation |
|--------|----------------------|---------------------|
| Timestep | Collapses to 1e-10 | Stable at ~1e-4 |
| Continuity error | >1e15 | <1e-9 |
| Completion | Diverges | Completes |

## Files

| File | Purpose |
|------|---------|
| `0/p` | Windkessel pressure BC (kinematic units, BDF3, implicit coupling) |
| `0/U` | Stabilised velocity BC (betaT=0.2, betaN=0.0) |
| `system/controlDict` | Solver control (loads libmodularWKPressure.so) |
| `system/fvSolution` | PIMPLE settings (pFinal=1.0 required) |
| `constant/boundaryData/inlet/` | Time-varying inlet velocity data |
