# OpenFOAM Windkessel (WK) Boundary Condition with Backflow Stabilization

A specialized OpenFOAM boundary condition implementing the three-element (RCR) Windkessel model for cardiovascular flow simulations with advanced backflow stabilization. This boundary condition is particularly useful for modeling the effects of downstream vasculature in cardiovascular CFD simulations, especially in complex 3D patient-specific geometries.

> **Current Status**: This implementation includes backflow stabilization based on Esmaily Moghadam et al. (2011), enabling stable simulations of complex arterial geometries with flow reversal. The stabilization has been tested on a coarctation of aorta (CoA) case, showing 14,000x improvement in timestep stability.

## Overview

The Windkessel boundary condition calculates outlet pressure based on flow rate using a lumped parameter model. This implementation provides:

- Three-element (RCR) Windkessel model support
- Multiple time discretization orders (1st, 2nd, and 3rd order)
- Robust state handling for case restarts
- Modular design for easy extension
- **NEW**: Backflow stabilization for complex 3D geometries
- **NEW**: Prevention of numerical instability at outlets with flow reversal

## Installation

### Prerequisites

- OpenFOAM v12.x
- C++ compiler compatible with your OpenFOAM installation
- wmake build system (included with OpenFOAM)

### Building

1. Source your OpenFOAM environment:
```bash
source /path/to/OpenFOAM-12/etc/bashrc
```

2. Navigate to the source directory and compile:
```bash
cd src/modularWKPressure
wmake
```

## Usage

### Boundary Condition Parameters

#### Pressure Boundary (modularWKPressure)

| Parameter | Description | Units |
|-----------|-------------|-------|
| `R` | Peripheral resistance | [m^-4] |
| `C` | Compliance | [m^4 s^2 kg^-1] |
| `Z` | Characteristic impedance | [m^-4] |
| `order` | Time discretization order (1-3) | - |
| `phi` | Name of the flux field (default: "phi") | - |
| `p0` | Initial pressure value | [m^2/s^2] |
| `value` | Initial uniform value | [m^2/s^2] |

#### Velocity Boundary (stabilizedWindkesselVelocity) - NEW

| Parameter | Description | Units | Default |
|-----------|-------------|-------|---------|
| `beta` | Stabilization coefficient | - | 1.0 |
| `enableStabilization` | Enable/disable backflow stabilization | - | true |

### Example Configuration

#### Basic Setup (without stabilization)

For simple 2D cases or when backflow is not expected, use only the pressure boundary in `0/p`:

```cpp
outlet
{
    type            modularWKPressure;
    phi             phi;
    order           2;
    R               1000;
    C               1e-6;
    Z               100;
    p0              0;
    value           uniform 0;
}
```

#### Stabilized Setup (recommended for 3D cases)

For complex 3D geometries with potential backflow, use both pressure and velocity boundaries:

In `0/p`:
```cpp
outlet
{
    type            modularWKPressure;
    phi             phi;
    order           2;
    R               1000;
    C               1e-6;
    Z               100;
    p0              0;
    value           uniform 0;
}
```

In `0/U`:
```cpp
outlet
{
    type                stabilizedWindkesselVelocity;
    beta                1.0;
    enableStabilization true;
}
```

### Initializing Historical Values

For a stable start of the simulation, it's important to properly initialize the historical values. The boundary condition maintains several historical values for both pressure and flow rate:

#### Required Values
- `p0`: Initial pressure value (required)
- `value`: Initial uniform pressure value (required)

#### Optional Historical Values
- `p_1`: Pressure at t-2dt (defaults to p0 if not specified)
- `q_1`: Flow rate at t-dt (required for all orders)
- `q_2`: Flow rate at t-2dt (required for 2nd and 3rd order, defaults to q_1)
- `q_3`: Flow rate at t-3dt (required for 3rd order, defaults to q_2)

Example configuration with all historical values for 3rd order accuracy:

```cpp
outlet
{
    type            modularWKPressure;
    phi             phi;
    order           3;
    R               1000;
    C               1e-6;
    Z               100;
    
    // Initial and historical pressures
    p0              0;      // Initial pressure
    p_1             0;      // Pressure at t-2dt
    value           uniform 0;
    
    // Historical flow rates
    q_1             0;      // Flow rate at t-dt
    q_2             0;      // Flow rate at t-2dt
    q_3             0;      // Flow rate at t-3dt
}
```

#### Initialization Strategies

1. **Cold Start** (when starting from rest):
   ```cpp
   outlet
   {
       type            modularWKPressure;
       phi             phi;
       order           1;      // Start with 1st order
       R               1000;
       C               1e-6;
       Z               100;
       p0              0;
       value           uniform 0;
       q_1             0;      // Zero initial flow
   }
   ```

2. **Steady State Start** (when starting from a known flow rate Q):
   ```cpp
   outlet
   {
       type            modularWKPressure;
       phi             phi;
       order           1;      // Start with 1st order
       R               1000;
       C               1e-6;
       Z               100;
       p0              #calc "$Q*$R";  // Initial pressure = Q*R
       value           uniform #calc "$Q*$R";
       q_1             $Q;     // Known flow rate
   }
   ```

3. **Gradual Order Increase** (recommended approach):
   - Start with 1st order for the first few time steps
   - Once stable, switch to 2nd order
   - Finally, switch to 3rd order if needed

Note: When restarting a simulation, the boundary condition automatically reads and uses the historical values from the previous run, ensuring continuity.

## Backflow Stabilization

### When to Use Stabilization

Backflow stabilization is crucial when:
- Running 3D patient-specific arterial geometries
- Outlets experience flow reversal during the cardiac cycle
- Simulation shows timestep collapse (deltaT < 1e-8)
- Continuity errors explode (> 1e10)

### Stabilization Parameters

- **beta**: Controls stabilization strength
  - 0.1-0.5: Conservative stabilization
  - 1.0: Recommended default
  - >1.0: Strong stabilization (use if still unstable)

### Troubleshooting Guide

| Symptom | Solution |
|---------|----------|
| Timestep < 1e-8 | Increase beta to 1.0-1.5 |
| Still unstable | Check mesh quality at outlets |
| Over-damped flow | Reduce beta to 0.5-0.8 |
| Slow convergence | Adjust solver tolerances |

## Tutorial Cases

### Basic 2D Case
The `tutorials/pitzDailyLESPulseWK/` demonstrates basic functionality using the standard OpenFOAM pitzDaily geometry.

### Complex 3D Case with Stabilization
The `tutorials/CoA_test/` demonstrates the stabilized boundary condition on a coarctation of aorta geometry with multiple outlets and complex flow patterns.

### Current Test Case
- Basic setup using pitzDaily geometry
- Integration with LES simulation
- Handling of pulsatile flow conditions

### Upcoming Validation Cases
- Patient-specific arterial geometries
- Physiologically relevant flow conditions
- Multiple outlet configurations
- Comparison with clinical data

To run the current tutorial:

```bash
cd tutorials/pitzDailyLESPulseWK
./Allrun
```

> **Note**: Stay tuned for updates as we add more physiologically relevant test cases and validation results.

## References

1. Westerhof, N., Lankhaar, J. W., & Westerhof, B. E. (2009). The arterial Windkessel. Medical & biological engineering & computing, 47(2), 131-141.
2. Esmaily Moghadam, M., Bazilevs, Y., Hsia, T. Y., Vignon-Clementel, I. E., & Marsden, A. L. (2011). A comparison of outlet boundary treatments for prevention of backflow divergence with relevance to blood flow simulations. Computational Mechanics, 48(3), 277-291.
3. OpenFOAM User Guide: [Boundary Conditions](https://www.openfoam.com/documentation/user-guide/4-boundaries/4.2-boundaries)

### Related Projects and Articles

This implementation was inspired by and builds upon previous work:

- [OpenFOAM-v8-Windkessel-code](https://github.com/EManchester/OpenFOAM-v8-Windkessel-code) by Emily Manchester et al. If you use this implementation or the original v8 version in your research, please cite:
  > Manchester, E. L., Pirola, S., Salmasi, M. Y., O'Regan, D. P., Athanasiou, T., and Xu, X. Y. (2021). Analysis of Turbulence Effects in a Patient-Specific Aorta with Aortic Valve Stenosis. Cardiovasc. Eng. Tech. 12, 438453. doi:10.1007/s13239-021-00536-9

- [Part 1: Modular Windkessel Boundary Condition in OpenFOAM v12](https://medium.com/@jiewang-share/part-1-modular-windkessel-boundary-condition-in-openfoam-v12-aaaab845923f) - A detailed tutorial on implementing and using the Windkessel boundary condition in OpenFOAM v12.

- [Part 2: Backflow Stabilization for Windkessel Boundary Conditions in OpenFOAM v12](https://medium.com/word-garden/part-2-backflow-stabilization-for-windkessel-boundary-conditions-in-openfoam-v12) - Advanced stabilization techniques for handling backflow in complex 3D cardiovascular simulations.
