# Pitz-Daily LES with Pulsatile Flow and Windkessel Outlet

This test case demonstrates the usage of the modular Windkessel boundary condition with a pulsatile flow in the Pitz-Daily backward-facing step geometry.

## Case Setup

- **Geometry**: Pitz-Daily backward-facing step
- **Flow**: Pulsatile inlet with LES turbulence modeling
- **Outlet**: 3-element Windkessel boundary condition
- **Time**: Transient simulation with second-order time discretization

## Directory Structure

```
pitzDailyLESPulseWK/
├── 0/                  # Initial conditions
│   ├── U              # Velocity
│   ├── p              # Pressure
│   ├── k              # Turbulent kinetic energy
│   ├── nut            # Turbulent viscosity
│   └── nuTilda        # Modified turbulent viscosity
├── constant/          # Physical properties
│   ├── transportProperties
│   └── turbulenceProperties
└── system/           # Solver settings
    ├── controlDict   # Time control and output
    ├── fvSchemes     # Discretization schemes
    └── fvSolution    # Linear solver settings
```

## Running the Case

1. Ensure the Windkessel boundary condition is compiled
2. Run the case:
   ```bash
   ./Allrun
   ```

## Post-processing

The case includes function objects for monitoring:
- Pressure and flow rate at the outlet
- Velocity profiles at various locations
- Turbulence statistics

Results can be visualized using:
```bash
paraFoam
``` 