# Modular Windkessel Pressure Boundary Condition with Backflow Stabilization

This is an implementation of a three-element Windkessel model as a pressure boundary condition for OpenFOAM, enhanced with backflow stabilization for complex 3D cardiovascular simulations. The boundary condition calculates outlet pressure based on flow rate using a lumped parameter model of the downstream vasculature, while the stabilization prevents numerical instability from flow reversal.

## Features

### Windkessel Pressure BC
- Three-element (RCR) Windkessel model
- Support for 1st, 2nd, and 3rd order time discretization
- Robust state handling for case restarts
- Modular design for easy extension

### Backflow Stabilization (NEW)
- Automatic backflow detection (v·n < 0)
- Esmaily Moghadam stabilization approach
- Prevents timestep collapse in complex geometries
- Tunable stabilization parameters

## Implementation Details

The boundary condition:
1. Reads the flow rate from the flux field
2. Solves the Windkessel ODE using the specified time discretization order
3. Updates the pressure boundary value
4. Stores historical values for the next timestep

### Parameters

#### Pressure BC (modularWKPressure)
- `R`: Peripheral resistance [m^-4]
- `C`: Compliance [m^4 s^2 kg^-1]
- `Z`: Characteristic impedance [m^-4]
- `order`: Time discretization order (1-3)
- `phi`: Name of the flux field (default: "phi")

#### Velocity BC (stabilizedWindkesselVelocity)
- `beta`: Stabilization coefficient (0.1-1.5, default: 1.0)
- `enableStabilization`: Enable/disable flag (default: true)

### State Variables

The boundary condition maintains historical values for:
- Pressure: p0_ (t-dt), p_1_ (t-2dt)
- Flow rate: q0_ (t-dt), q_1_ (t-2dt), q_2_ (t-3dt), q_3_ (t-4dt)

## Compilation

```bash
wmake
```

## Usage Examples

### Basic Setup (2D cases)

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

### Stabilized Setup (3D cases with backflow)

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

## Performance Impact

Example results from coarctation of aorta case:
- **Without stabilization**: deltaT = 1e-10 s (collapsed)
- **With stabilization**: deltaT = 1e-4 s (stable)
- **Improvement factor**: 14,000x 