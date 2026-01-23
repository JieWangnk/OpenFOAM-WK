# modularWKPressure Library

**Advanced outlet boundary conditions for cardiovascular CFD simulations in OpenFOAM 12**

This library provides three boundary condition types for physiologically-accurate outlet modeling:

1. **modularWKPressure** - Three-element Windkessel (RCR) pressure BC
2. **vectorFittingImpedance** - Multi-pole rational function impedance BC
3. **stabilizedWindkesselVelocity** - Backflow-stabilized velocity BC

---

## Installation

```bash
cd ~/OpenFOAM/<user>-12/src/modularWKPressure
wmake libso
```

Verify installation:
```bash
ls $FOAM_USER_LIBBIN/libmodularWKPressure.so
```

Load in `system/controlDict`:
```cpp
libs ("libmodularWKPressure.so");
```

---

## Boundary Conditions

### 1. modularWKPressure

Three-element Windkessel model for outlet pressure using **kinematic units**.

**Circuit analogy:**
```
          Z           R
Inlet ----/\/\/\-----/\/\/\----- Ground (p=0)
                  |
                 === C
                  |
                 --- Ground
```

**Parameters:**
| Parameter | Unit | Description |
|-----------|------|-------------|
| R | m⁻¹·s⁻¹ | Peripheral resistance (kinematic) |
| C | m·s² | Compliance (kinematic) |
| Z | m⁻¹·s⁻¹ | Characteristic impedance (kinematic) |
| p0 | m²/s² | Reference pressure (kinematic) |
| order | - | Time discretization (1, 2, or 3) |
| couplingMode | - | `explicit` or `implicit` (recommended) |

**Unit conversion (ρ = 1060 kg/m³):**
- R_kin = R_dyn / ρ  (Pa·s/m³ → m⁻¹·s⁻¹)
- C_kin = C_dyn × ρ  (m³/Pa → m·s²)
- p_kin = p_dyn / ρ  (Pa → m²/s²)

**Example (`0/p`):**
```cpp
outlet1
{
    type            modularWKPressure;
    phi             phi;
    U               U;
    couplingMode    implicit;
    order           3;

    // Kinematic units
    R               1108115.88;     // m⁻¹·s⁻¹
    C               9.024e-07;      // m·s²
    Z               268937.17;      // m⁻¹·s⁻¹
    p0              10.06;          // m²/s² (~80 mmHg)

    // State variables (for restart)
    p_1             10.06;
    q_1             0;
    q_2             0;
    q_3             0;

    value           uniform 10.06;
}
```

### 2. vectorFittingImpedance

Multi-pole rational function impedance model using recursive convolution.

**Mathematical model:**
```
Z(s) = d + Σᵢ rᵢ/(s - pᵢ)
```

**Parameters:**
| Parameter | Unit | Description |
|-----------|------|-------------|
| nPoles | - | Number of poles (typically 4-6). `order` accepted for backward compat. |
| directTerm | Pa·s/m³ or m⁻¹·s⁻¹ | Direct feedthrough term |
| poles | rad/s | Pole locations (must be negative) |
| residues | Pa/m³ or m⁻¹ | Pole residues |
| rho | kg/m³ | Fluid density (for unit conversion) |
| impedanceUnits | - | `dynamic` (default) or `kinematic` |
| couplingMode | - | `explicit` or `implicit` |

**Example (`0/p`):**
```cpp
outlet1
{
    type                vectorFittingImpedance;
    phi                 phi;
    U                   U;
    couplingMode        implicit;
    nPoles              4;

    // Dynamic units (Pa-based) - internally converted using rho
    impedanceUnits      dynamic;
    rho                 1060.0;

    directTerm          6.93e+05;
    poles               (-6.71 -12.24 -53.35 -76.91);
    residues            (-4.80e+06 1.46e+06 1.23e+08 -2.15e+08);

    value               uniform 10.06;
}
```

To use kinematic units directly (no conversion):
```cpp
    impedanceUnits      kinematic;
    // Parameters already in m⁻¹·s⁻¹ and m⁻¹
```

### 3. stabilizedWindkesselVelocity

Backflow-stabilized velocity BC using directional (`directionMixed`) control.

**Two-parameter control:**
- **betaT** - Tangential damping (vortex suppression)
- **betaN** - Normal damping (flow reversal control)

**Physics:** During backflow (phi < 0):
- Tangential velocity is damped toward zero (betaT controls strength)
- Normal velocity responds to pressure (betaN=0 for Windkessel compatibility)

**Formula:**
```
valueFraction = H(-phi) × [betaN × (n⊗n) + betaT × (I - n⊗n)]
```

**Parameters:**
| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| betaT | 0-1 | 0.2 | Tangential stabilization |
| betaN | 0-1 | 0.0 | Normal stabilization |
| beta | 0-1 | - | Legacy (maps to betaT) |
| enableStabilization | bool | true | Enable/disable |
| dampingFactor | 0-1 | 1.0 | Additional multiplier |

**Recommended values:**
| Scenario | betaT | betaN |
|----------|-------|-------|
| Mild backflow | 0.1-0.2 | 0.0 |
| Moderate backflow | 0.2-0.3 | 0.0 |
| Severe backflow | 0.3-0.5 | 0.0-0.1 |

**Example (`0/U`):**
```cpp
outlet1
{
    type                  stabilizedWindkesselVelocity;
    phi                   phi;
    betaT                 0.2;    // Tangential damping
    betaN                 0.0;    // Normal: free for Windkessel
    enableStabilization   true;
    value                 uniform (0 0 0);
}
```

---

## Complete Outlet Setup

Pair pressure and velocity BCs on each outlet:

**`0/p`:**
```cpp
outlet1
{
    type            modularWKPressure;
    phi             phi;
    U               U;
    couplingMode    implicit;
    order           3;
    R               1108115.88;
    C               9.024e-07;
    Z               268937.17;
    p0              10.06;
    value           uniform 10.06;
}
```

**`0/U`:**
```cpp
outlet1
{
    type                  stabilizedWindkesselVelocity;
    phi                   phi;
    betaT                 0.2;
    betaN                 0.0;
    value                 uniform (0 0 0);
}
```

---

## Typical Pressure Ranges

| Pressure | Dynamic (Pa) | Kinematic (m²/s²) |
|----------|-------------|-------------------|
| 80 mmHg (diastolic) | 10,666 | 10.06 |
| 120 mmHg (systolic) | 16,000 | 15.09 |
| 100 mmHg (mean) | 13,333 | 12.58 |

---

## Parallel Execution

All BCs use `gSum()` for global flux summation. Safe for any decomposition method.

---

## Restart Behavior

State variables are written to time directories:
- modularWKPressure: `p0_`, `p_1_`, `p_2_`, `q_1_`, `q_2_`, `q_3_`
- vectorFittingImpedance: `stateVariables_`, `q_1_`

**Clean restart after re-decomposition:**
1. `reconstructPar`
2. Delete `processor*` directories
3. `decomposePar`

---

## Utilities

Python tools for impedance extraction and vector fitting:

| File | Purpose |
|------|---------|
| `impedanceVectorFit.py` | Fit rational function to impedance data |
| `extractImpedanceFromCFD.py` | Extract Z(f) from CFD results |
| `convert_postProcessing_to_impedance.py` | Process OpenFOAM postProcessing output |

See `utilities/README_VECTOR_FITTING.md` for workflow details.

---

## Troubleshooting

**Pressure values ~10,000 instead of ~10:**
Using dynamic units (Pa) instead of kinematic. Convert: p_kin = p_dyn / rho

**Simulation diverges during diastole:**
Increase betaT (0.2 → 0.3 → 0.5)

**Timestep collapses:**
- Reduce betaT if too aggressive
- Check mesh quality at outlets
- Verify PIMPLE convergence

**Missing phi field error:**
Ensure using incompressible solver (pimpleFoam, simpleFoam) that creates phi

---

## References

1. Westerhof et al. (2009). "The arterial Windkessel." Med Biol Eng Comput, 47:131-141.
2. Esmaily Moghadam et al. (2011). "A comparison of outlet boundary treatments." Ann Biomed Eng, 39(5):1384-1399.
3. Gustavsen & Semlyen (1999). "Rational approximation of frequency domain responses." IEEE Trans Power Deliv, 14(3):1052-1061.

---

## License

GNU General Public License v3.0 (same as OpenFOAM)
