#!/usr/bin/env python3
"""
Extract Impedance from OpenFOAM CFD Results
============================================

Post-processes OpenFOAM simulation results to compute frequency-domain
impedance Z(ω) = P̂(ω)/Q̂(ω) for use with vector fitting.

Workflow:
    1. Extract time-series pressure P(t) and flowrate Q(t) from outlets
    2. Apply FFT to get frequency-domain P̂(ω) and Q̂(ω)
    3. Compute impedance Z(ω) = P̂(ω)/Q̂(ω)
    4. Save to CSV for input to impedanceVectorFit.py

Usage:
    # Extract from OpenFOAM case directory:
    python extractImpedanceFromCFD.py -case /path/to/openfoam_case -outlet outlet1

    # Process multiple outlets:
    python extractImpedanceFromCFD.py -case ./openfoam_implicit -outlet outlet1 outlet2

    # With custom time range (exclude transients):
    python extractImpedanceFromCFD.py -case ./case -outlet outlet1 -start-time 0.5

Author: OpenFOAM-WK Development Team
Date: 2025-11-15
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import sys
import re


class OpenFOAMCaseReader:
    """Read OpenFOAM case data for impedance extraction."""

    def __init__(self, case_dir):
        """
        Initialize reader for OpenFOAM case.

        Parameters
        ----------
        case_dir : str or Path
            Path to OpenFOAM case directory
        """
        self.case_dir = Path(case_dir)
        if not self.case_dir.exists():
            raise ValueError(f"Case directory not found: {case_dir}")

        # Detect time directories
        self.time_dirs = self._get_time_directories()
        print(f"Found {len(self.time_dirs)} time directories")

    def _get_time_directories(self):
        """Get sorted list of time directories."""
        time_dirs = []
        for item in self.case_dir.iterdir():
            if item.is_dir():
                try:
                    time_value = float(item.name)
                    time_dirs.append((time_value, item))
                except ValueError:
                    # Not a numeric directory (e.g., "0.orig", "constant", "system")
                    pass

        # Sort by time value
        time_dirs.sort(key=lambda x: x[0])
        return time_dirs

    def extract_outlet_timeseries(self, outlet_name, start_time=0.0, rho=1060.0):
        """
        Extract pressure and flowrate time series for an outlet.

        Parameters
        ----------
        outlet_name : str
            Name of the outlet patch
        start_time : float
            Start time (to exclude initial transients)
        rho : float
            Fluid density [kg/m³] to convert kinematic → dynamic pressure

        Returns
        -------
        times : ndarray
            Time points [s]
        pressures : ndarray
            Mean outlet pressure [Pa] (dynamic)
        flowrates : ndarray
            Outlet flowrate [m³/s]
        """
        times = []
        pressures = []
        flowrates = []

        for time_value, time_dir in self.time_dirs:
            if time_value < start_time:
                continue

            # Read pressure field (kinematic [m²/s²])
            p_file = time_dir / 'p'
            if not p_file.exists():
                print(f"Warning: p file not found in {time_dir}")
                continue

            p_patch = self._read_patch_field(p_file, outlet_name)
            if p_patch is not None:
                # Convert kinematic → dynamic pressure
                p_mean = np.mean(p_patch) * rho  # [Pa]
                pressures.append(p_mean)
            else:
                print(f"Warning: patch {outlet_name} not found in {p_file}")
                continue

            # Read flux field phi [m³/s]
            phi_file = time_dir / 'phi'
            if not phi_file.exists():
                print(f"Warning: phi file not found in {time_dir}")
                continue

            phi_patch = self._read_patch_field(phi_file, outlet_name)
            if phi_patch is not None:
                q = np.sum(phi_patch)  # Total flowrate [m³/s]
                flowrates.append(q)
            else:
                print(f"Warning: patch {outlet_name} not found in {phi_file}")
                continue

            times.append(time_value)

        if len(times) == 0:
            raise ValueError(f"No data extracted for outlet {outlet_name}")

        return np.array(times), np.array(pressures), np.array(flowrates)

    def _read_patch_field(self, field_file, patch_name):
        """
        Read OpenFOAM patch field data.

        This is a simplified parser - for complex cases use PyFoam or similar.
        """
        with open(field_file, 'r') as f:
            content = f.read()

        # Find the boundaryField section
        boundary_match = re.search(
            r'boundaryField\s*\{(.*?)\n\}',
            content,
            re.DOTALL
        )

        if not boundary_match:
            return None

        boundary_content = boundary_match.group(1)

        # Find the specific patch
        patch_pattern = rf'{patch_name}\s*\{{(.*?)\n\s*\}}'
        patch_match = re.search(patch_pattern, boundary_content, re.DOTALL)

        if not patch_match:
            return None

        patch_content = patch_match.group(1)

        # Extract uniform or nonuniform values
        # Uniform case: "value uniform 10.06;"
        uniform_match = re.search(r'value\s+uniform\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)',
                                  patch_content)
        if uniform_match:
            value = float(uniform_match.group(1))
            return np.array([value])

        # Nonuniform case: "value nonuniform List<scalar> n(val1 val2 ...);"
        nonuniform_match = re.search(
            r'value\s+nonuniform\s+List<\w+>\s*\n?\s*\d+\s*\(([\s\S]*?)\)',
            patch_content
        )
        if nonuniform_match:
            values_str = nonuniform_match.group(1)
            values = np.array([float(x) for x in values_str.split()])
            return values

        return None


def compute_impedance_fft(times, pressures, flowrates, window='hann'):
    """
    Compute frequency-domain impedance using FFT.

    Z(ω) = P̂(ω) / Q̂(ω)

    Parameters
    ----------
    times : ndarray
        Time points [s]
    pressures : ndarray
        Pressure time series [Pa]
    flowrates : ndarray
        Flowrate time series [m³/s]
    window : str
        Window function ('hann', 'hamming', 'blackman', or None)

    Returns
    -------
    frequencies : ndarray
        Frequency points [Hz] (positive frequencies only)
    impedance_complex : ndarray
        Complex impedance [Pa·s/m³]
    """
    # Check for uniform sampling
    dt = np.diff(times)
    if np.std(dt) / np.mean(dt) > 0.01:
        print("Warning: Non-uniform time sampling detected. Interpolating...")
        # Interpolate to uniform grid
        t_uniform = np.linspace(times[0], times[-1], len(times))
        pressures = np.interp(t_uniform, times, pressures)
        flowrates = np.interp(t_uniform, times, flowrates)
        times = t_uniform
        dt = np.mean(np.diff(times))
    else:
        dt = np.mean(dt)

    # Apply window to reduce spectral leakage
    if window:
        if window == 'hann':
            win = np.hanning(len(times))
        elif window == 'hamming':
            win = np.hamming(len(times))
        elif window == 'blackman':
            win = np.blackman(len(times))
        else:
            raise ValueError(f"Unknown window: {window}")

        pressures = pressures * win
        flowrates = flowrates * win

    # Compute FFT
    P_fft = np.fft.fft(pressures)
    Q_fft = np.fft.fft(flowrates)

    # Compute impedance Z = P/Q
    # Handle division by zero
    Z_fft = np.zeros_like(P_fft, dtype=complex)
    valid = np.abs(Q_fft) > 1e-15
    Z_fft[valid] = P_fft[valid] / Q_fft[valid]

    # Get positive frequencies only
    n = len(times)
    freqs = np.fft.fftfreq(n, dt)
    positive_freq_idx = freqs > 0

    frequencies = freqs[positive_freq_idx]
    impedance_complex = Z_fft[positive_freq_idx]

    return frequencies, impedance_complex


def save_impedance_csv(filename, frequencies, impedance_complex):
    """
    Save impedance data to CSV for input to impedanceVectorFit.py.

    Parameters
    ----------
    filename : str
        Output CSV file path
    frequencies : ndarray
        Frequency points [Hz]
    impedance_complex : ndarray
        Complex impedance [Pa·s/m³]
    """
    df = pd.DataFrame({
        'frequency': frequencies,
        'magnitude': np.abs(impedance_complex),
        'phase': np.angle(impedance_complex),
        'real': impedance_complex.real,
        'imag': impedance_complex.imag
    })

    df.to_csv(filename, index=False)
    print(f"Impedance data saved to {filename}")
    print(f"  {len(frequencies)} frequency points")
    print(f"  Range: {frequencies[0]:.3f} - {frequencies[-1]:.3f} Hz")


def plot_timeseries(times, pressures, flowrates, save_path=None):
    """Plot pressure and flowrate time series."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.plot(times, pressures / 133.322, linewidth=2)  # Convert to mmHg
    ax1.set_xlabel('Time [s]')
    ax1.set_ylabel('Pressure [mmHg]')
    ax1.set_title('Outlet Pressure Time Series')
    ax1.grid(True, alpha=0.3)

    ax2.plot(times, flowrates * 1e6, linewidth=2)  # Convert to mL/s
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Flowrate [mL/s]')
    ax2.set_title('Outlet Flowrate Time Series')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Time series plot saved to {save_path}")
    else:
        plt.show()


def plot_impedance(frequencies, impedance_complex, save_path=None):
    """Plot impedance Bode diagrams."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Magnitude
    ax1.loglog(frequencies, np.abs(impedance_complex), linewidth=2)
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('Impedance Magnitude [Pa·s/m³]')
    ax1.set_title('Impedance Spectrum: Magnitude')
    ax1.grid(True, which='both', alpha=0.3)

    # Phase
    ax2.semilogx(frequencies, np.rad2deg(np.angle(impedance_complex)), linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Phase [degrees]')
    ax2.set_title('Impedance Spectrum: Phase')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Impedance plot saved to {save_path}")
    else:
        plt.show()


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Extract impedance from OpenFOAM CFD results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract impedance for outlet1:
  python extractImpedanceFromCFD.py -case ./openfoam_implicit -outlet outlet1

  # Multiple outlets:
  python extractImpedanceFromCFD.py -case ./case -outlet outlet1 outlet2 outlet3

  # Exclude initial transients (start from t=0.5s):
  python extractImpedanceFromCFD.py -case ./case -outlet outlet1 -start-time 0.5

Output:
  - impedance_<outlet>.csv: Impedance data for vector fitting
  - timeseries_<outlet>.png: Pressure and flowrate plots
  - impedance_<outlet>.png: Bode diagrams
        """
    )

    parser.add_argument('-case', '-c', required=True,
                        help='OpenFOAM case directory')
    parser.add_argument('-outlet', '-o', nargs='+', required=True,
                        help='Outlet patch name(s)')
    parser.add_argument('-start-time', type=float, default=0.0,
                        help='Start time to exclude transients [s], default=0')
    parser.add_argument('-rho', type=float, default=1060.0,
                        help='Fluid density [kg/m³], default=1060 (blood)')
    parser.add_argument('-window', choices=['hann', 'hamming', 'blackman', 'none'],
                        default='hann',
                        help='FFT window function, default=hann')
    parser.add_argument('-no-plots', action='store_true',
                        help='Skip plot generation')

    args = parser.parse_args()

    window = args.window if args.window != 'none' else None

    print(f"\n{'='*60}")
    print(f"Extract Impedance from OpenFOAM CFD Results")
    print(f"{'='*60}\n")

    # Initialize reader
    try:
        reader = OpenFOAMCaseReader(args.case)
    except Exception as e:
        print(f"Error reading case: {e}")
        sys.exit(1)

    # Process each outlet
    for outlet_name in args.outlet:
        print(f"\n{'='*60}")
        print(f"Processing outlet: {outlet_name}")
        print(f"{'='*60}\n")

        try:
            # Extract time series
            times, pressures, flowrates = reader.extract_outlet_timeseries(
                outlet_name,
                start_time=args.start_time,
                rho=args.rho
            )

            print(f"Extracted {len(times)} time points")
            print(f"  Time range: {times[0]:.3f} - {times[-1]:.3f} s")
            print(f"  Duration: {times[-1] - times[0]:.3f} s")
            print(f"  Mean pressure: {np.mean(pressures):.2f} Pa "
                  f"({np.mean(pressures)/133.322:.2f} mmHg)")
            print(f"  Mean flowrate: {np.mean(flowrates)*1e6:.4f} mL/s")

            # Compute impedance
            frequencies, impedance_complex = compute_impedance_fft(
                times, pressures, flowrates, window=window
            )

            print(f"\nComputed impedance spectrum:")
            print(f"  {len(frequencies)} frequency points")
            print(f"  Fundamental frequency: {frequencies[0]:.3f} Hz")
            print(f"  Max frequency: {frequencies[-1]:.3f} Hz")

            # Save CSV
            csv_filename = f'impedance_{outlet_name}.csv'
            save_impedance_csv(csv_filename, frequencies, impedance_complex)

            # Generate plots
            if not args.no_plots:
                plot_timeseries(times, pressures, flowrates,
                               save_path=f'timeseries_{outlet_name}.png')
                plot_impedance(frequencies, impedance_complex,
                              save_path=f'impedance_{outlet_name}.png')

            print(f"\n✓ Successfully processed {outlet_name}")
            print(f"  Next step: Run vector fitting with:")
            print(f"    python impedanceVectorFit.py -input {csv_filename} -order 4\n")

        except Exception as e:
            print(f"Error processing {outlet_name}: {e}")
            import traceback
            traceback.print_exc()
            continue

    print(f"\n{'='*60}")
    print(f"Extraction complete!")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
