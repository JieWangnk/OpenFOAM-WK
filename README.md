# OpenFOAM Windkessel (WK) Boundary Condition

A specialized OpenFOAM boundary condition implementing the three-element (RCR) Windkessel model for cardiovascular flow simulations. This boundary condition is particularly useful for modeling the effects of downstream vasculature in cardiovascular CFD simulations.

## Overview

The Windkessel boundary condition calculates outlet pressure based on flow rate using a lumped parameter model. This implementation provides:

- Three-element (RCR) Windkessel model support
- Multiple time discretization orders (1st, 2nd, and 3rd order)
- Robust state handling for case restarts
- Modular design for easy extension

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

- `R`: Peripheral resistance [m^-4]
- `C`: Compliance [m^4 s^2 kg^-1]
- `Z`: Characteristic impedance [m^-4]
- `order`: Time discretization order (1-3)
- `phi`: Name of the flux field (default: "phi")
- `p0`: Initial pressure value
- `value`: Initial uniform value

### Example Configuration

Add the following to your case's boundary field (e.g., `0/p`):

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

## Tutorial Case

A tutorial case demonstrating the usage of the Windkessel boundary condition is provided in `tutorials/pitzDailyLESPulseWK/`. This case showcases:

- Basic setup of the Windkessel boundary condition
- Integration with LES simulation
- Handling of pulsatile flow conditions

To run the tutorial:

```bash
cd tutorials/pitzDailyLESPulseWK
./Allrun
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit pull requests or open issues for:

- Bug fixes
- Feature enhancements
- Documentation improvements
- Additional tutorial cases

## References

1. Westerhof, N., Lankhaar, J. W., & Westerhof, B. E. (2009). The arterial Windkessel. Medical & biological engineering & computing, 47(2), 131-141.
2. OpenFOAM User Guide: [Boundary Conditions](https://www.openfoam.com/documentation/user-guide/4-boundaries/4.2-boundaries) 