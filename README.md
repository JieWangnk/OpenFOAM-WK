# OpenFOAM Windkessel (WK) Boundary Condition with Backflow Stabilization

A specialized OpenFOAM boundary condition implementing the three-element (RCR) Windkessel model for cardiovascular flow simulations with advanced two-parameter backflow stabilization. This boundary condition is particularly useful for modeling the effects of downstream vasculature in cardiovascular CFD simulations, especially in complex 3D patient-specific geometries.

> **Current Status**: This implementation includes two-parameter directional backflow stabilization (βT, βN) based on Esmaily Moghadam et al. (2011), enabling stable simulations of complex arterial geometries with flow reversal while preserving Windkessel pressure-flow coupling.

## Publication

A paper describing this software is in preparation for **SoftwareX** (Elsevier):

> **modularWKPressure: A Two-Parameter Backflow-Stabilised Windkessel Boundary Condition for Cardiovascular CFD in OpenFOAM**

The paper includes:
- Verification tests (ODE accuracy, CFD coupling, impedance response)
- Patient-specific aortic application comparing stabilisation strategies
- Formal pass/fail criteria for haemodynamic neutrality

Paper draft and supplementary materials: [`OpenFOAM_WK_paper/`](../OpenFOAM_WK_paper/)

---

## Three-Element Windkessel Model

**Windkessel ODE:**

$$P = ZQ + P_c, \quad \frac{dP_c}{dt} = \frac{Q}{C} - \frac{P_c}{RC}$$

**Unit conversions (kinematic ↔ dynamic):**

| Parameter | Kinematic (OpenFOAM) | Dynamic (Physical) | Conversion |
|-----------|---------------------|-------------------|------------|
| Pressure p | m²·s⁻² | Pa | p_kin = p_dyn / ρ |
| Resistance R, Z | m⁻¹·s⁻¹ | Pa·s·m⁻³ | R_kin = R_dyn / ρ |
| Compliance C | m·s² | m³·Pa⁻¹ | C_kin = C_dyn · ρ |

With ρ = 1060 kg/m³ for blood.

---

## Overview

The Windkessel boundary condition calculates outlet pressure based on flow rate using a lumped parameter model. This implementation provides:

- Three-element (RCR) Windkessel model support
- Multiple time discretization orders (BDF1, BDF2, BDF3)
- Robust state handling for case restarts
- Modular design for easy extension
- **Two-parameter backflow stabilization** (βT for tangential, βN for normal)
- Prevention of numerical instability at outlets with flow reversal

---

## Installation

### Prerequisites

- OpenFOAM v12.x
- C++ compiler compatible with your OpenFOAM installation
- wmake build system (included with OpenFOAM)

### Building

```bash
# Source OpenFOAM environment
source /opt/openfoam12/etc/bashrc

# Build the library (portable path)
cd $WM_PROJECT_USER_DIR/src/modularWKPressure
wmake libso

# Verify library exists
ls $FOAM_USER_LIBBIN/libmodularWKPressure.so
```

---

## Usage

### Boundary Condition Parameters

#### Pressure Boundary (`modularWKPressure`)

| Parameter | Description | Units (kinematic) |
|-----------|-------------|-------------------|
| `R` | Peripheral resistance | m⁻¹·s⁻¹ |
| `C` | Compliance | m·s² |
| `Z` | Characteristic impedance | m⁻¹·s⁻¹ |
| `order` | BDF time discretization order (1-3) | - |
| `phi` | Name of flux field | - |
| `p0` | Initial capacitor pressure | m²/s² |
| `value` | Initial uniform value | m²/s² |

#### Velocity Boundary (`stabilizedWindkesselVelocity`)

| Parameter | Description | Range | Default |
|-----------|-------------|-------|---------|
| `betaT` | Tangential damping coefficient | 0.0–1.0 | 0.3 |
| `betaN` | Normal damping coefficient | 0.0–1.0 | 0.0 |
| `phi` | Name of flux field | - | phi |

**Recommended settings:**
- **βT = 0.3**: Suppresses tangential vortices during backflow
- **βN = 0.0**: Preserves Windkessel pressure-flow coupling (keep at zero for RCR outlets)

### Example Configuration

#### Recommended Setup (two-parameter stabilization)

In `0/p`:
```cpp
outlet
{
    type            modularWKPressure;
    phi             phi;
    order           3;          // BDF3 for accuracy
    R               1.11e6;     // [m⁻¹·s⁻¹] kinematic
    C               9.02e-7;    // [m·s²] kinematic
    Z               2.69e5;     // [m⁻¹·s⁻¹] kinematic
    p0              0;
    value           uniform 0;
}
```

In `0/U`:
```cpp
outlet
{
    type            stabilizedWindkesselVelocity;
    phi             phi;
    betaT           0.3;        // Tangential damping (recommended)
    betaN           0.0;        // Normal: free for Windkessel coupling
    value           uniform (0 0 0);
}
```

In `system/controlDict`:
```cpp
libs ( "libmodularWKPressure.so" );
```

#### Basic Setup (without stabilization)

For simple 2D cases or when backflow is not expected:

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
    type            zeroGradient;
}
```

---

## Two-Parameter Backflow Stabilization

### Tensor Formulation

$$\mathbf{F} = H(-\phi) \left( \beta_N \mathbf{n}\mathbf{n}^T + \beta_T (\mathbf{I} - \mathbf{n}\mathbf{n}^T) \right)$$

where:
- **n** is the outward unit normal
- **I** is the identity tensor
- **H(−φ)** is the Heaviside backflow indicator

### Why βN = 0 by Default?

**βN = 0 is the physically correct default for Windkessel coupling.**

| Direction | Parameter | Default | Rationale |
|-----------|-----------|---------|-----------|
| Tangential | βT | 0.2–0.3 | Suppresses vortices; doesn't affect P-Q relationship |
| **Normal** | βN | **0.0** | Must be free to respond to Windkessel pressure |

Setting βN > 0 constrains normal velocity toward zero during backflow, which **fights against** the pressure model trying to drive reverse flow.

### When to Use βN > 0

βN is a **safety valve**, not the default:

1. Truncated outlets with strong numerical reflections
2. Severe backflow that βT alone cannot stabilise
3. Under-resolved outlet regions (coarse mesh, poor orthogonality)

**Typical range if needed:** βN = 0.02–0.10

### Troubleshooting Guide

| Symptom | Solution |
|---------|----------|
| Timestep collapse (< 1e-8) | Increase βT to 0.3–0.5 |
| Still unstable | Add small βN (0.05) as last resort |
| Over-damped flow | Reduce βT to 0.1–0.2 |
| Slow convergence | Check mesh quality at outlets |

---

## Verification

### Verification Metrics

**Backflow volume fraction:**

$$\eta_{bf}(t) = \frac{|Q^{-}(t)|}{Q^{+}(t) + |Q^{-}(t)|}$$

**Haemodynamic neutrality:**

$$\Delta\bar{Q} = \frac{\bar{Q}_{stab} - \bar{Q}_{ref}}{\bar{Q}_{ref}}, \quad \Delta\bar{p} = \frac{\bar{p}_{stab} - \bar{p}_{ref}}{\bar{p}_{ref}}$$

### Pass/Fail Criteria

| Test | Metric | Pass Criterion |
|------|--------|----------------|
| RCR ODE accuracy | BDF3 vs analytical | Error < 5% |
| Steady-state coupling | P vs (Z+R)×Q | Error < 5% |
| Impedance response | \|Z(ω)\| at f_heart | Error < 5%, phase < 5° |
| No-harm condition | Stabilised vs reference (forward flow) | \|Δp̄\| < 1%, \|ΔQ̄\| < 1% |
| βN sensitivity | βN=0.05 vs βN=0.0 | \|Δp̄\| < 2%, \|ΔQ̄\| < 2% |

---

## Initializing Historical Values

The boundary condition maintains historical values for BDF schemes:

| Value | Description | Required for |
|-------|-------------|--------------|
| `p0` | Initial capacitor pressure | All orders |
| `p_1` | Pressure at t-2Δt | Order ≥ 2 |
| `q_1` | Flow rate at t-Δt | All orders |
| `q_2` | Flow rate at t-2Δt | Order ≥ 2 |
| `q_3` | Flow rate at t-3Δt | Order 3 |

Example for BDF3:
```cpp
outlet
{
    type            modularWKPressure;
    phi             phi;
    order           3;
    R               1.11e6;
    C               9.02e-7;
    Z               2.69e5;

    p0              0;
    p_1             0;
    q_1             0;
    q_2             0;
    q_3             0;
    value           uniform 0;
}
```

---

## Tutorial Cases

### Basic 2D Case
`tutorials/pitzDailyLESPulseWK/` - Basic functionality using pitzDaily geometry.

### Complex 3D Case
`tutorials/CoA_test/` - Coarctation of aorta with multiple outlets and stabilization.

```bash
cd tutorials/pitzDailyLESPulseWK
./Allrun
```

---

## Methods Description (Paper-Ready)

### Outlet boundary conditions: RCR pressure model with backflow-stabilised velocity

Outlet haemodynamics were represented using a three-element Windkessel (RCR) pressure boundary condition coupled with a backflow-stabilised velocity treatment. For incompressible simulations, OpenFOAM's pressure variable is the **kinematic pressure** *p = p_dyn / ρ* (units m²·s⁻²). The outlet flow rate was computed from the face flux field *φ* as the patch integral

$$Q(t) = \int_{\Gamma_{\text{out}}} \mathbf{U} \cdot \mathbf{n} \, dA \approx \sum_{f \in \Gamma_{\text{out}}} \phi_f$$

using a global reduction in parallel runs to ensure processor-independent *Q*.

#### Three-element Windkessel (RCR) pressure outlet

At each outlet, the Windkessel model was written in dynamic form as *P = ZQ + P_c*, where *Z* is the proximal resistance and *P_c* is the capacitor pressure satisfying

$$\frac{dP_c}{dt} = \frac{Q}{C} - \frac{P_c}{RC}$$

with distal resistance *R* and compliance *C*. The boundary pressure applied to the solver was the kinematic form *p = P/ρ*, with kinematic parameters defined by *R_kin = R_dyn/ρ*, *Z_kin = Z_dyn/ρ*, and *C_kin = C_dyn·ρ*. The Windkessel ODE was advanced once per CFD time step using a backward differentiation scheme (BDF1–BDF3), and the resulting outlet pressure was held fixed during PIMPLE outer iterations to improve robustness for stiff outlet parameters and pulsatile conditions.

#### Backflow stabilisation at outlets

To suppress non-physical re-entry at truncated outlets and improve convergence, a directional mixed velocity condition was applied only during backflow. On each outlet face, a tensor weight was defined as

$$\mathbf{F} = H(-\phi) \left( \beta_N \mathbf{n}\mathbf{n}^T + \beta_T (\mathbf{I} - \mathbf{n}\mathbf{n}^T) \right)$$

where **n** is the outward unit normal, **I** is the identity tensor, and *H(−φ)* is a Heaviside indicator based on the sign of the local face flux *φ* (with a small deadband around *φ = 0* to avoid switching due to numerical noise). The tangential coefficient *β_T* damps tangential velocity components during reverse flow, while *β_N* controls damping in the normal direction; for Windkessel outlets *β_N* was set to zero by default to avoid constraining the normal outflow that is determined by the pressure–flow model. The stabilisation was implemented through the `directionMixed` formulation so that the boundary contribution enters the momentum system consistently during linearisation, and no explicit patch overwrite (`evaluate()`) was performed during coefficient updates.

---

## Conclusion

> **The recommended configuration (βT = 0.3, βN = 0) reduces backflow indicators and improves solver robustness while preserving mean outlet pressure and flow within predefined tolerances (|Δp̄| < 1%, |ΔQ̄| < 1%). The βN parameter is kept at zero to maintain consistency with the Windkessel pressure-flow relationship; sensitivity tests confirm that small normal damping (βN = 0.05–0.10) does not materially alter integral outlet haemodynamics (< 2% change).**

---

## Citation

If you use this software in your research, please cite:

```bibtex
@article{wang2026modularwkpressure,
  title={modularWKPressure: A Two-Parameter Backflow-Stabilised Windkessel
         Boundary Condition for Cardiovascular CFD in OpenFOAM},
  author={Wang, Jie and others},
  journal={SoftwareX},
  year={2026},
  publisher={Elsevier}
}
```

---

## References

1. Westerhof, N., Lankhaar, J. W., & Westerhof, B. E. (2009). The arterial Windkessel. Medical & biological engineering & computing, 47(2), 131-141.
2. Esmaily Moghadam, M., Bazilevs, Y., Hsia, T. Y., Vignon-Clementel, I. E., & Marsden, A. L. (2011). A comparison of outlet boundary treatments for prevention of backflow divergence with relevance to blood flow simulations. Computational Mechanics, 48(3), 277-291.
3. OpenFOAM User Guide: [Boundary Conditions](https://www.openfoam.com/documentation/user-guide/4-boundaries/4.2-boundaries)

### Related Projects and Articles

This implementation was inspired by and builds upon previous work:

- [OpenFOAM-v8-Windkessel-code](https://github.com/EManchester/OpenFOAM-v8-Windkessel-code) by Emily Manchester et al. If you use this implementation or the original v8 version in your research, please cite:
  > Manchester, E. L., Pirola, S., Salmasi, M. Y., O'Regan, D. P., Athanasiou, T., and Xu, X. Y. (2021). Analysis of Turbulence Effects in a Patient-Specific Aorta with Aortic Valve Stenosis. Cardiovasc. Eng. Tech. 12, 438–453. doi:10.1007/s13239-021-00536-9

- [Part 1: Modular Windkessel Boundary Condition in OpenFOAM v12](https://medium.com/@jiewang-share/part-1-modular-windkessel-boundary-condition-in-openfoam-v12-aaaab845923f) - A detailed tutorial on implementing and using the Windkessel boundary condition in OpenFOAM v12.

- [Part 2: Backflow Stabilization for Windkessel Boundary Conditions in OpenFOAM v12](https://medium.com/word-garden/part-2-backflow-stabilization-for-windkessel-boundary-conditions-in-openfoam-v12) - Advanced stabilization techniques for handling backflow in complex 3D cardiovascular simulations.

---

**January 2026** | OpenFOAM 12 | modularWKPressure
