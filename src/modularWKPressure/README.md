# Modular Windkessel Pressure Boundary Condition

This is an implementation of a three-element Windkessel model as a pressure boundary condition for OpenFOAM. The boundary condition calculates outlet pressure based on flow rate using a lumped parameter model of the downstream vasculature.

## Features

- Three-element (RCR) Windkessel model
- Support for 1st, 2nd, and 3rd order time discretization
- Robust state handling for case restarts
- Modular design for easy extension

## Implementation Details

The boundary condition:
1. Reads the flow rate from the flux field
2. Solves the Windkessel ODE using the specified time discretization order
3. Updates the pressure boundary value
4. Stores historical values for the next timestep

### Parameters

- `R`: Peripheral resistance [m^-4]
- `C`: Compliance [m^4 s^2 kg^-1]
- `Z`: Characteristic impedance [m^-4]
- `order`: Time discretization order (1-3)
- `phi`: Name of the flux field (default: "phi")

### State Variables

The boundary condition maintains historical values for:
- Pressure: p0_ (t-dt), p_1_ (t-2dt)
- Flow rate: q0_ (t-dt), q_1_ (t-2dt), q_2_ (t-3dt), q_3_ (t-4dt)

## Compilation

```bash
wmake
```

## Usage Example

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