#!/usr/bin/env python3
"""
Convert OpenFOAM postProcessing Data to Impedance CSV
======================================================

This script reads .dat files from postProcessing/surfaceFieldValue/
and converts them to impedance CSV files compatible with impedanceVectorFit.py

Usage:
    python convert_postProcessing_to_impedance.py \
        -case /path/to/case \
        -outlet outlet1 outlet2 outlet3 \
        -start-time 0.5 \
        -rho 1060

Input:
    postProcessing/{outlet}Pressure/0/surfaceFieldValue.dat
    postProcessing/{outlet}Flow/0/surfaceFieldValue.dat

Output:
    impedance_{outlet}.csv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
from pathlib import Path


def read_surface_field_value(dat_file):
    """
    Read OpenFOAM surfaceFieldValue .dat file

    Format:
    # Time    average(p) or sum(phi)
    0.0      10665.76
    0.05     10723.45
    ...

    Returns:
        time (array), value (array)
    """
    try:
        data = np.loadtxt(dat_file, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        time = data[:, 0]
        value = data[:, 1]
        return time, value
    except Exception as e:
        raise RuntimeError(f"Failed to read {dat_file}: {e}")


def compute_impedance_fft(time, pressure, flow_rate, start_time=0.0, rho=1060.0):
    """
    Compute impedance from pressure and flow rate time series using FFT

    Args:
        time: Time array [s]
        pressure: Pressure array [Pa]
        flow_rate: Flow rate array [m³/s]
        start_time: Start time for analysis (exclude transients)
        rho: Fluid density [kg/m³]

    Returns:
        freq: Frequency array [Hz]
        Z_real: Real part of impedance [Pa·s/m³]
        Z_imag: Imaginary part of impedance [Pa·s/m³]
    """
    # Trim transient period
    mask = time >= start_time
    time = time[mask]
    pressure = pressure[mask]
    flow_rate = flow_rate[mask]

    if len(time) < 10:
        raise ValueError(f"Not enough data points after start_time={start_time}s (only {len(time)} points)")

    # Check sampling
    dt = np.diff(time)
    dt_mean = np.mean(dt)
    dt_std = np.std(dt)

    print(f"  Time points: {len(time)}")
    print(f"  Time range: {time[0]:.3f} - {time[-1]:.3f} s")
    print(f"  Mean dt: {dt_mean:.6f} s (std: {dt_std:.6f} s)")

    # Check if uniformly sampled
    if dt_std / dt_mean > 0.1:
        print(f"  Warning: Non-uniform time sampling detected!")
        print(f"           FFT results may be inaccurate")
        print(f"           Consider interpolating to uniform grid")

    # Compute FFT
    N = len(time)
    duration = time[-1] - time[0]

    # Apply window to reduce spectral leakage
    window = np.hanning(N)
    pressure_windowed = (pressure - np.mean(pressure)) * window
    flow_windowed = (flow_rate - np.mean(flow_rate)) * window

    # FFT
    P_fft = np.fft.rfft(pressure_windowed)
    Q_fft = np.fft.rfft(flow_windowed)
    freq = np.fft.rfftfreq(N, dt_mean)

    # Impedance Z = P / Q (in frequency domain)
    # Add small epsilon to avoid division by zero
    epsilon = 1e-10 * np.max(np.abs(Q_fft))
    Z_complex = P_fft / (Q_fft + epsilon)

    # Remove DC component (freq = 0)
    freq = freq[1:]
    Z_complex = Z_complex[1:]

    # Only keep physiologically relevant frequencies (0.1 - 20 Hz)
    mask = (freq >= 0.1) & (freq <= 20.0)
    freq = freq[mask]
    Z_complex = Z_complex[mask]

    print(f"  Frequency range: {freq[0]:.3f} - {freq[-1]:.3f} Hz")
    print(f"  Number of frequencies: {len(freq)}")

    return freq, Z_complex.real, Z_complex.imag


def save_impedance_csv(filename, freq, Z_real, Z_imag):
    """
    Save impedance to CSV file in format expected by impedanceVectorFit.py
    """
    Z_mag = np.sqrt(Z_real**2 + Z_imag**2)
    Z_phase = np.degrees(np.arctan2(Z_imag, Z_real))

    df = pd.DataFrame({
        'frequency': freq,  # Note: impedanceVectorFit.py expects 'frequency' or 'freq'
        'real': Z_real,
        'imag': Z_imag,
        'magnitude': Z_mag,
        'phase': Z_phase
    })

    df.to_csv(filename, index=False)
    print(f"  Saved: {filename}")


def plot_time_series(time, pressure, flow_rate, outlet_name):
    """
    Plot pressure and flow rate time series
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

    # Pressure
    ax1.plot(time, pressure / 133.322, 'b-', linewidth=1)
    ax1.set_ylabel('Pressure [mmHg]')
    ax1.set_title(f'{outlet_name} - Time Series')
    ax1.grid(True, alpha=0.3)

    # Flow rate
    ax2.plot(time, flow_rate * 1e6, 'r-', linewidth=1)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Flow Rate [mL/s]')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    filename = f'timeseries_{outlet_name}.png'
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"  Saved: {filename}")


def plot_impedance(freq, Z_real, Z_imag, outlet_name):
    """
    Plot impedance magnitude and phase
    """
    Z_mag = np.sqrt(Z_real**2 + Z_imag**2)
    Z_phase = np.degrees(np.arctan2(Z_imag, Z_real))

    # Convert to CGS units for plotting (dyn·s/cm⁵)
    Z_mag_cgs = Z_mag * 0.1  # Pa·s/m³ -> dyn·s/cm⁵

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6))

    # Magnitude
    ax1.loglog(freq, Z_mag_cgs, 'b-', linewidth=2, label='Extracted')
    ax1.set_ylabel('|Z| [dyn·s/cm⁵]')
    ax1.set_title(f'{outlet_name} - Impedance Spectrum')
    ax1.grid(True, which='both', alpha=0.3)
    ax1.legend()

    # Phase
    ax2.semilogx(freq, Z_phase, 'r-', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Phase [degrees]')
    ax2.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    filename = f'impedance_{outlet_name}.png'
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"  Saved: {filename}")


def process_outlet(case_dir, outlet_name, start_time, rho, plot=True):
    """
    Process one outlet: read postProcessing data, compute impedance, save CSV
    """
    print(f"\nProcessing {outlet_name}...")

    # Find postProcessing directories
    postproc_dir = Path(case_dir) / 'postProcessing'

    # Look for pressure data
    pressure_dirs = list(postproc_dir.glob(f'{outlet_name}Pressure/*/surfaceFieldValue.dat'))
    if not pressure_dirs:
        raise FileNotFoundError(f"No pressure data found for {outlet_name} in {postproc_dir}")

    pressure_file = pressure_dirs[0]
    print(f"  Reading pressure: {pressure_file}")

    # Look for flow data
    flow_dirs = list(postproc_dir.glob(f'{outlet_name}Flow/*/surfaceFieldValue.dat'))
    if not flow_dirs:
        raise FileNotFoundError(f"No flow data found for {outlet_name} in {postproc_dir}")

    flow_file = flow_dirs[0]
    print(f"  Reading flow: {flow_file}")

    # Read data
    time_p, pressure = read_surface_field_value(pressure_file)
    time_q, flow_rate = read_surface_field_value(flow_file)

    # Check if time arrays match (check length first to avoid numpy broadcast error)
    arrays_match = (len(time_p) == len(time_q)) and np.allclose(time_p, time_q, atol=1e-6)

    if not arrays_match:
        print(f"  Warning: Pressure and flow time arrays don't match")
        print(f"           Pressure: {len(time_p)} points from {time_p[0]:.3f} to {time_p[-1]:.3f} s")
        print(f"           Flow:     {len(time_q)} points from {time_q[0]:.3f} to {time_q[-1]:.3f} s")

        # Find overlapping time range
        t_min = max(time_p[0], time_q[0])
        t_max = min(time_p[-1], time_q[-1])

        print(f"           Using overlapping range: {t_min:.3f} - {t_max:.3f} s")

        # Use the longer time array as base and interpolate the shorter one
        if len(time_p) >= len(time_q):
            time = time_p[(time_p >= t_min) & (time_p <= t_max)]
            pressure = pressure[(time_p >= t_min) & (time_p <= t_max)]
            # Interpolate flow to match
            flow_rate = np.interp(time, time_q, flow_rate)
        else:
            time = time_q[(time_q >= t_min) & (time_q <= t_max)]
            flow_rate = flow_rate[(time_q >= t_min) & (time_q <= t_max)]
            # Interpolate pressure to match
            pressure = np.interp(time, time_p, pressure)

        print(f"           Aligned to {len(time)} points")
    else:
        time = time_p

    # Plot time series
    if plot:
        plot_time_series(time, pressure, flow_rate, outlet_name)

    # Compute impedance
    print(f"  Computing impedance via FFT...")
    freq, Z_real, Z_imag = compute_impedance_fft(
        time, pressure, flow_rate,
        start_time=start_time,
        rho=rho
    )

    # Save CSV
    csv_filename = f'impedance_{outlet_name}.csv'
    save_impedance_csv(csv_filename, freq, Z_real, Z_imag)

    # Plot impedance
    if plot:
        plot_impedance(freq, Z_real, Z_imag, outlet_name)

    print(f"  ✓ {outlet_name} complete")

    return freq, Z_real, Z_imag


def main():
    parser = argparse.ArgumentParser(
        description='Convert OpenFOAM postProcessing data to impedance CSV',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Process all outlets
    python convert_postProcessing_to_impedance.py \\
        -case ~/GitHub/AortaCFD-app/output/BPM120/run_20251115_131246/openfoam_stable \\
        -outlet outlet1 outlet2 outlet3 outlet4 \\
        -start-time 0.5 \\
        -rho 1060

    # Process single outlet without plots
    python convert_postProcessing_to_impedance.py \\
        -case /path/to/case \\
        -outlet outlet1 \\
        -no-plots
"""
    )

    parser.add_argument('-case', required=True,
                       help='Path to OpenFOAM case directory')
    parser.add_argument('-outlet', nargs='+', required=True,
                       help='Outlet patch names (e.g., outlet1 outlet2)')
    parser.add_argument('-start-time', type=float, default=0.5,
                       help='Start time for analysis (exclude transients) [default: 0.5s]')
    parser.add_argument('-rho', type=float, default=1060.0,
                       help='Fluid density [kg/m³] [default: 1060]')
    parser.add_argument('-no-plots', action='store_true',
                       help='Skip plotting (faster)')

    args = parser.parse_args()

    # Validate case directory
    case_dir = Path(args.case)
    if not case_dir.exists():
        print(f"Error: Case directory not found: {case_dir}")
        return 1

    postproc_dir = case_dir / 'postProcessing'
    if not postproc_dir.exists():
        print(f"Error: postProcessing directory not found: {postproc_dir}")
        print(f"")
        print(f"You need to run extract_using_openfoam.sh first:")
        print(f"  bash extract_using_openfoam.sh {case_dir} {args.start_time}")
        return 1

    print("=" * 72)
    print("Convert OpenFOAM postProcessing Data to Impedance CSV")
    print("=" * 72)
    print(f"")
    print(f"Case:       {case_dir}")
    print(f"Outlets:    {', '.join(args.outlet)}")
    print(f"Start time: {args.start_time} s")
    print(f"Density:    {args.rho} kg/m³")
    print(f"")

    # Process each outlet
    for outlet in args.outlet:
        try:
            process_outlet(
                case_dir=case_dir,
                outlet_name=outlet,
                start_time=args.start_time,
                rho=args.rho,
                plot=not args.no_plots
            )
        except Exception as e:
            import traceback
            print(f"  ✗ Error processing {outlet}: {e}")
            if args.start_time < 0:  # Debug mode trigger (use negative start_time to enable)
                print("\nFull traceback:")
                traceback.print_exc()
            continue

    print("")
    print("=" * 72)
    print("Conversion complete!")
    print("=" * 72)
    print("")
    print("Generated files:")
    print("  - impedance_*.csv (for vector fitting)")
    if not args.no_plots:
        print("  - timeseries_*.png (pressure and flow vs time)")
        print("  - impedance_*.png (impedance spectrum)")
    print("")
    print("Next step: Run vector fitting")
    print("  python impedanceVectorFit.py -input impedance_outlet1.csv -order 4 ...")
    print("")

    return 0


if __name__ == '__main__':
    exit(main())
