#!/usr/bin/env python3
"""
Vector Fitting for Impedance Boundary Conditions
=================================================

Estimates vectorFittingImpedance BC parameters from frequency-domain impedance data
using the Vector Fitting algorithm (Gustavsen & Semlyen, 1999).

Mathematical Model:
    Z(s) = d + Σᵢ₌₁ᴺ rᵢ/(s - pᵢ)

where:
    s     = Laplace variable (s = iω)
    d     = direct feedthrough term [Pa·s/m³]
    rᵢ    = residues [Pa/m³]  (NOT Pa·s/m³!)
    pᵢ    = poles [rad/s], must be negative for stability
    N     = order (typically 4-6, recommended 4)

Input Data Formats:
    1. CSV with columns: frequency [Hz], magnitude [Pa·s/m³], phase [rad]
    2. CSV with columns: frequency [Hz], real [Pa·s/m³], imag [Pa·s/m³]
    3. Clinical impedance spectra from 4D Flow MRI
    4. 1D model impedance output

Output:
    - OpenFOAM dictionary entry for vectorFittingImpedance BC
    - Validation plots comparing fitted vs. original impedance
    - Pole-residue diagram
    - Error metrics (RMSE, max error)

Usage:
    # From frequency-magnitude-phase CSV:
    python impedanceVectorFit.py -input data/impedance_outlet1.csv -order 4

    # With custom output:
    python impedanceVectorFit.py -input data/impedance.csv -order 6 -output myBC

    # From clinical data:
    python impedanceVectorFit.py -input clinical/4DFlow_impedance.mat -order 4

References:
    - Gustavsen & Semlyen (1999), "Rational approximation of frequency
      domain responses by vector fitting," IEEE Trans. Power Delivery
    - Fevola et al. (2023), "A vector fitting approach for the automated
      estimation of lumped boundary conditions of 1D circulation models,"
      Frontiers in Physiology, 14:1250204

Author: OpenFOAM-WK Development Team
Date: 2025-11-15
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal
from scipy.optimize import least_squares
import argparse
import sys
from pathlib import Path


class VectorFittingImpedance:
    """
    Vector Fitting algorithm for impedance boundary conditions.

    Fits a rational function Z(s) = d + Σrᵢ/(s-pᵢ) to frequency-domain
    impedance data using an iterative pole relocation algorithm.
    """

    def __init__(self, order=4, max_iterations=20, tolerance=1e-4):
        """
        Initialize Vector Fitting solver.

        Parameters
        ----------
        order : int
            Number of poles (typically 4-6, recommended 4)
        max_iterations : int
            Maximum iterations for pole relocation
        tolerance : float
            Convergence tolerance for relative error
        """
        self.order = order
        self.max_iterations = max_iterations
        self.tolerance = tolerance

        # Results (to be filled by fit())
        self.poles = None
        self.residues = None
        self.direct_term = None
        self.fitted = False

    def fit(self, frequencies, impedance_complex):
        """
        Fit rational function to impedance data using Vector Fitting.

        Parameters
        ----------
        frequencies : array_like
            Frequency points [Hz]
        impedance_complex : array_like
            Complex impedance Z(f) = Z_real + j·Z_imag [Pa·s/m³]

        Returns
        -------
        dict
            Fitted parameters: poles, residues, direct_term
        """
        # Convert frequency to angular frequency
        omega = 2 * np.pi * frequencies  # [rad/s]
        s = 1j * omega  # Laplace variable

        # Initialize poles (logarithmically spaced negative real poles)
        omega_min = omega[0] if omega[0] > 0 else omega[1]
        omega_max = omega[-1]
        initial_poles = -np.logspace(
            np.log10(omega_min),
            np.log10(omega_max),
            self.order
        )

        print(f"Vector Fitting: order={self.order}, n_freq={len(frequencies)}")
        print(f"Initial poles: {initial_poles}")

        # Iterative pole relocation
        poles = initial_poles.copy()
        for iteration in range(self.max_iterations):
            # Step 1: Solve for residues given current poles
            residues, direct_term = self._solve_residues(s, impedance_complex, poles)

            # Step 2: Pole relocation (simplified - use stability-preserving update)
            poles_new = self._relocate_poles(s, impedance_complex, poles, residues, direct_term)

            # Check convergence
            pole_change = np.max(np.abs(poles_new - poles) / np.abs(poles))
            if pole_change < self.tolerance:
                print(f"Converged in {iteration+1} iterations (pole change: {pole_change:.2e})")
                break

            poles = poles_new

            if iteration == self.max_iterations - 1:
                print(f"Warning: Max iterations reached. Final pole change: {pole_change:.2e}")

        # Store results
        self.poles = poles
        self.residues = residues
        self.direct_term = direct_term
        self.fitted = True

        # Validate stability
        self._validate_stability()

        return {
            'poles': self.poles,
            'residues': self.residues,
            'direct_term': self.direct_term
        }

    def _solve_residues(self, s, Z_data, poles):
        """
        Solve for residues given poles using linear least squares.

        Z(s) = d + Σ rᵢ/(s - pᵢ)

        This is linear in [d, r₁, r₂, ..., rₙ].
        """
        n_freq = len(s)
        n_poles = len(poles)

        # Build design matrix A where Z ≈ A·x
        # x = [d, r₁, r₂, ..., rₙ]ᵀ
        A = np.zeros((n_freq, n_poles + 1), dtype=complex)
        A[:, 0] = 1.0  # Direct term column

        for i, pole in enumerate(poles):
            A[:, i+1] = 1.0 / (s - pole)  # Residue columns

        # Solve: min ||Z_data - A·x||²
        # Use real formulation for numerical stability
        A_real = np.vstack([A.real, A.imag])
        Z_real = np.hstack([Z_data.real, Z_data.imag])

        x, residuals, rank, singular = np.linalg.lstsq(A_real, Z_real, rcond=None)

        direct_term = x[0]
        residues = x[1:]

        return residues, direct_term

    def _relocate_poles(self, s, Z_data, poles, residues, direct_term):
        """
        Relocate poles for next iteration (simplified version).

        Full Vector Fitting uses a sophisticated pole relocation scheme.
        This simplified version uses a gradient-based update.
        """
        # For simplicity, keep poles on negative real axis and adjust magnitudes
        # Full implementation would solve a separate fitting problem

        # Compute current fit error using current poles and residues
        Z_fit = np.zeros(len(s), dtype=complex) + direct_term
        for pole, residue in zip(poles, residues):
            Z_fit += residue / (s - pole)
        error = np.abs(Z_data - Z_fit)

        # Adjust poles based on error distribution
        # Poles with high local error should move to those frequency regions
        weights = np.zeros(len(poles))
        for i in range(len(poles)):
            # Weight by proximity to high-error regions
            omega = np.abs(s.imag)
            pole_omega = np.abs(poles[i])
            distance = np.abs(omega - pole_omega)
            weights[i] = np.sum(error / (distance + 1e-10))

        # Adjust pole positions (keep on negative real axis)
        adjustment_factor = 1.1
        poles_new = poles * (1.0 + 0.1 * weights / np.max(weights + 1e-10))

        # Ensure stability (all poles negative)
        poles_new = -np.abs(poles_new)

        return poles_new

    def evaluate(self, frequencies, with_imag=False):
        """
        Evaluate fitted impedance at given frequencies.

        Parameters
        ----------
        frequencies : array_like
            Frequency points [Hz]
        with_imag : bool
            If True, return complex impedance. If False, return magnitude.

        Returns
        -------
        Z : array_like
            Impedance (complex if with_imag=True, magnitude otherwise)
        """
        if not self.fitted:
            raise RuntimeError("Must call fit() before evaluate()")

        omega = 2 * np.pi * frequencies
        s = 1j * omega

        Z = self.direct_term * np.ones_like(s)
        for pole, residue in zip(self.poles, self.residues):
            Z += residue / (s - pole)

        if with_imag:
            return Z
        else:
            return np.abs(Z)

    def _validate_stability(self):
        """Validate that all poles are negative (stability requirement)."""
        if np.any(self.poles >= 0):
            unstable_indices = np.where(self.poles >= 0)[0]
            print(f"WARNING: Unstable poles detected at indices {unstable_indices}")
            print(f"Pole values: {self.poles[unstable_indices]}")
            print("Force-correcting to negative values...")
            self.poles[unstable_indices] = -np.abs(self.poles[unstable_indices])

    def plot_fit(self, frequencies, impedance_complex, save_path=None):
        """
        Plot fitted vs. original impedance (Bode plots).

        Parameters
        ----------
        frequencies : array_like
            Original frequency points [Hz]
        impedance_complex : array_like
            Original complex impedance [Pa·s/m³]
        save_path : str, optional
            Path to save plot (if None, display only)
        """
        Z_fitted = self.evaluate(frequencies, with_imag=True)

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Magnitude plot
        ax1.loglog(frequencies, np.abs(impedance_complex), 'o',
                   label='Original', markersize=6, alpha=0.7)
        ax1.loglog(frequencies, np.abs(Z_fitted), '-',
                   label=f'Fitted (order {self.order})', linewidth=2)
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Impedance Magnitude [Pa·s/m³]')
        ax1.set_title('Impedance Bode Plot: Magnitude')
        ax1.grid(True, which='both', alpha=0.3)
        ax1.legend()

        # Phase plot
        phase_orig = np.angle(impedance_complex)
        phase_fit = np.angle(Z_fitted)
        ax2.semilogx(frequencies, np.rad2deg(phase_orig), 'o',
                     label='Original', markersize=6, alpha=0.7)
        ax2.semilogx(frequencies, np.rad2deg(phase_fit), '-',
                     label=f'Fitted (order {self.order})', linewidth=2)
        ax2.set_xlabel('Frequency [Hz]')
        ax2.set_ylabel('Phase [degrees]')
        ax2.set_title('Impedance Bode Plot: Phase')
        ax2.grid(True, which='both', alpha=0.3)
        ax2.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Bode plot saved to {save_path}")
        else:
            plt.show()

    def plot_poles(self, save_path=None):
        """
        Plot pole-residue diagram.

        Parameters
        ----------
        save_path : str, optional
            Path to save plot
        """
        fig, ax = plt.subplots(figsize=(8, 6))

        # Plot poles in complex plane
        ax.plot(self.poles.real, self.poles.imag, 'x',
                markersize=12, markeredgewidth=2, label='Poles')

        # Mark residue magnitudes
        for i, (pole, res) in enumerate(zip(self.poles, self.residues)):
            ax.annotate(f'r{i+1}={res:.2e}',
                       xy=(pole.real, pole.imag),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8)

        ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
        ax.axvline(x=0, color='r', linestyle='--', alpha=0.5, label='Stability boundary')
        ax.set_xlabel('Real(s) [rad/s]')
        ax.set_ylabel('Imag(s) [rad/s]')
        ax.set_title('Pole-Residue Diagram')
        ax.grid(True, alpha=0.3)
        ax.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Pole diagram saved to {save_path}")
        else:
            plt.show()

    def write_openfoam_bc(self, outlet_name='outlet', rho=1060.0, p0=13332.0):
        """
        Generate OpenFOAM boundary condition dictionary entry.

        Parameters
        ----------
        outlet_name : str
            Name of the outlet patch
        rho : float
            Fluid density [kg/m³], default 1060 for blood
        p0 : float
            Reference pressure [Pa], default 13332 (100 mmHg)

        Returns
        -------
        str
            OpenFOAM dictionary entry
        """
        if not self.fitted:
            raise RuntimeError("Must call fit() before write_openfoam_bc()")

        # Convert kinematic initial pressure
        p0_kin = p0 / rho

        bc_text = f"""    {outlet_name}
    {{
        type                  vectorFittingImpedance;
        phi                   phi;
        U                     U;
        couplingMode          implicit;  // Use implicit for stability
        order                 {self.order};

        // Vector fitting parameters (DYNAMIC units - auto-converted)
        directTerm            {self.direct_term:.10e};  // [Pa·s/m³]
        poles                 ({' '.join(f'{p:.10e}' for p in self.poles)});  // [rad/s]
        residues              ({' '.join(f'{r:.10e}' for r in self.residues)});  // [Pa/m³]

        // Fluid properties
        rho                   {rho};  // [kg/m³]

        // Initial conditions
        value                 uniform {p0_kin:.6f};  // [m²/s²] kinematic = {p0}/rho
    }}
"""
        return bc_text


def load_impedance_data(file_path):
    """
    Load impedance data from CSV file.

    Supports multiple formats:
    1. frequency, magnitude, phase
    2. frequency, real, imag
    3. frequency, impedance (assumes real)

    Parameters
    ----------
    file_path : str
        Path to CSV file

    Returns
    -------
    frequencies : ndarray
        Frequency points [Hz]
    impedance_complex : ndarray
        Complex impedance [Pa·s/m³]
    """
    df = pd.read_csv(file_path)

    # Detect format
    columns = [c.lower() for c in df.columns]

    if 'frequency' not in columns and 'freq' not in columns:
        raise ValueError("CSV must have 'frequency' or 'freq' column")

    freq_col = 'frequency' if 'frequency' in columns else 'freq'
    frequencies = df[freq_col].values

    if 'magnitude' in columns and 'phase' in columns:
        # Format 1: magnitude-phase
        magnitude = df['magnitude'].values
        phase = df['phase'].values  # radians
        impedance_complex = magnitude * np.exp(1j * phase)
        print(f"Loaded impedance: magnitude-phase format ({len(frequencies)} points)")

    elif 'real' in columns and 'imag' in columns:
        # Format 2: real-imaginary
        real = df['real'].values
        imag = df['imag'].values
        impedance_complex = real + 1j * imag
        print(f"Loaded impedance: real-imag format ({len(frequencies)} points)")

    elif 'impedance' in columns:
        # Format 3: real-valued impedance only
        impedance_complex = df['impedance'].values + 0j
        print(f"Loaded impedance: real-only format ({len(frequencies)} points)")
        print("Warning: Phase information not available, assuming zero phase")

    else:
        raise ValueError(
            "CSV must have either:\n"
            "  - 'magnitude' and 'phase' columns, or\n"
            "  - 'real' and 'imag' columns, or\n"
            "  - 'impedance' column"
        )

    return frequencies, impedance_complex


def main():
    """Main entry point for command-line usage."""
    parser = argparse.ArgumentParser(
        description='Vector Fitting for Impedance Boundary Conditions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Fit 4th-order model from CSV data:
  python impedanceVectorFit.py -input data/impedance_outlet1.csv -order 4

  # Fit 6th-order model with custom output:
  python impedanceVectorFit.py -input impedance.csv -order 6 -output outlet2_BC

  # Generate plots only:
  python impedanceVectorFit.py -input impedance.csv -order 4 -plots-only

Input CSV format (choose one):
  1. frequency[Hz], magnitude[Pa·s/m³], phase[rad]
  2. frequency[Hz], real[Pa·s/m³], imag[Pa·s/m³]
  3. frequency[Hz], impedance[Pa·s/m³]  (real-valued)
        """
    )

    parser.add_argument('-input', '-i', required=True,
                        help='Input CSV file with impedance data')
    parser.add_argument('-order', '-o', type=int, default=4,
                        help='Vector fitting order (number of poles), default=4')
    parser.add_argument('-output', default='vectorFitting_BC.txt',
                        help='Output file for OpenFOAM BC, default=vectorFitting_BC.txt')
    parser.add_argument('-rho', type=float, default=1060.0,
                        help='Fluid density [kg/m³], default=1060 (blood)')
    parser.add_argument('-p0', type=float, default=13332.0,
                        help='Reference pressure [Pa], default=13332 (100 mmHg)')
    parser.add_argument('-outlet-name', default='outlet',
                        help='Outlet patch name, default=outlet')
    parser.add_argument('-plots-only', action='store_true',
                        help='Generate plots without writing BC file')
    parser.add_argument('-no-plots', action='store_true',
                        help='Skip plot generation')

    args = parser.parse_args()

    # Validate inputs
    if not Path(args.input).exists():
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)

    if args.order < 1 or args.order > 10:
        print(f"Warning: Unusual order={args.order}. Recommended range: 4-6")

    # Load data
    print(f"\n{'='*60}")
    print(f"Vector Fitting for Impedance Boundary Conditions")
    print(f"{'='*60}\n")

    try:
        frequencies, impedance_complex = load_impedance_data(args.input)
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)

    # Print data summary
    print(f"\nData summary:")
    print(f"  Frequency range: {frequencies[0]:.3f} - {frequencies[-1]:.3f} Hz")
    print(f"  Number of points: {len(frequencies)}")
    print(f"  Impedance range: {np.min(np.abs(impedance_complex)):.2e} - "
          f"{np.max(np.abs(impedance_complex)):.2e} Pa·s/m³")

    # Perform vector fitting
    print(f"\n{'='*60}")
    print(f"Performing Vector Fitting (order={args.order})...")
    print(f"{'='*60}\n")

    vf = VectorFittingImpedance(order=args.order)
    params = vf.fit(frequencies, impedance_complex)

    # Display results
    print(f"\n{'='*60}")
    print(f"Vector Fitting Results")
    print(f"{'='*60}\n")
    print(f"Direct term d = {params['direct_term']:.6e} Pa·s/m³\n")
    print(f"Poles [rad/s]:")
    for i, p in enumerate(params['poles']):
        print(f"  p{i+1} = {p:.6e}")
    print(f"\nResidues [Pa/m³]:")
    for i, r in enumerate(params['residues']):
        print(f"  r{i+1} = {r:.6e}")

    # Compute fit quality
    Z_fitted = vf.evaluate(frequencies, with_imag=True)
    rmse = np.sqrt(np.mean(np.abs(impedance_complex - Z_fitted)**2))
    max_error = np.max(np.abs(impedance_complex - Z_fitted))
    rel_rmse = rmse / np.mean(np.abs(impedance_complex)) * 100

    print(f"\nFit quality:")
    print(f"  RMSE: {rmse:.4e} Pa·s/m³ ({rel_rmse:.2f}%)")
    print(f"  Max error: {max_error:.4e} Pa·s/m³")

    # Generate plots
    if not args.no_plots:
        print(f"\n{'='*60}")
        print(f"Generating plots...")
        print(f"{'='*60}\n")

        output_stem = Path(args.output).stem
        vf.plot_fit(frequencies, impedance_complex,
                    save_path=f'{output_stem}_bode.png')
        vf.plot_poles(save_path=f'{output_stem}_poles.png')

    # Write OpenFOAM BC
    if not args.plots_only:
        bc_text = vf.write_openfoam_bc(
            outlet_name=args.outlet_name,
            rho=args.rho,
            p0=args.p0
        )

        with open(args.output, 'w') as f:
            f.write(bc_text)

        print(f"\n{'='*60}")
        print(f"OpenFOAM BC written to: {args.output}")
        print(f"{'='*60}\n")
        print(bc_text)

    print(f"\nVector fitting complete!")


if __name__ == '__main__':
    main()
