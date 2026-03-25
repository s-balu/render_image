#!/usr/bin/env python3
"""
get_snap_dt.py — compute the snap_dt parameter for render_image velocity interpolation.

snap_dt is the time interval between two snapshots in units of
(simulation length unit) / (simulation velocity unit).

For Gadget-2/4, AREPO, and SWIFT HDF5 snapshots:
  - Positions are in comoving Mpc/h
  - Velocities are peculiar velocities with a sqrt(a) factor:
      v_stored = sqrt(a) * dx/dt   [km/s]
    The physical peculiar velocity is v_pec = v_stored / sqrt(a).

  - To convert stored velocity to comoving displacement per unit snap_dt:
      dx_comoving = v_stored * snap_dt / a
    where snap_dt has units of (Mpc/h) / (km/s).

Usage
-----
  python3 get_snap_dt.py --snap_a1 0.5 --snap_a2 0.52 \\
      --H0 67.74 --Omega_m 0.309 --Omega_L 0.691

  python3 get_snap_dt.py --snap_z1 1.0 --snap_z2 0.9 \\
      --H0 67.74 --Omega_m 0.309 --Omega_L 0.691 --h 0.6774

  python3 get_snap_dt.py --hdf5 snapshot_050.hdf5 --hdf5_next snapshot_051.hdf5

If two HDF5 files are provided, the cosmological parameters and scale
factors are read directly from the file headers.
"""

import argparse
import math
import sys

# ── Physical constants ─────────────────────────────────────────────────────────
# 1 Mpc = 3.085677581e19 km
# H0 in units of km/s/Mpc → H0 * (km/s/Mpc) * (Mpc / 3.086e19 km) = H0/3.086e19 s^-1
KM_PER_MPC = 3.085677581e19   # km per Mpc


def H_over_H0(a, Omega_m, Omega_L, Omega_k=0.0):
    """E(a) = H(a)/H0 for flat or curved ΛCDM."""
    return math.sqrt(Omega_m / a**3 + Omega_k / a**2 + Omega_L)


def comoving_time_integral(a1, a2, H0, Omega_m, Omega_L, n_steps=10000):
    """
    Compute ∫_{a1}^{a2} da / (a * H(a))  in seconds.

    This is the physical time elapsed between scale factors a1 and a2.
    H(a) = H0 * E(a),  H0 in km/s/Mpc.
    """
    # Convert H0 to s^-1
    H0_si = H0 / KM_PER_MPC   # s^-1

    da = (a2 - a1) / n_steps
    total = 0.0
    for i in range(n_steps):
        a = a1 + (i + 0.5) * da
        total += 1.0 / (a * H_over_H0(a, Omega_m, Omega_L) * H0_si)
    return total * da   # seconds


def snap_dt_gadget(a1, a2, H0, Omega_m, Omega_L, h_little, n_steps=10000):
    """
    Compute snap_dt in units of (Mpc/h) / (km/s) for Gadget-convention
    velocity interpolation.

    Gadget stores velocities as  v_stored = sqrt(a) * a * (dx_comoving/dt)
    where x_comoving is in Mpc/h and t is cosmic time in seconds.

    To shift comoving position by one snapshot interval:
        dx_comoving = v_stored / (sqrt(a) * a) * dt_physical
                    = v_stored * snap_dt / a           (our convention)
    where snap_dt = dt_physical [s] * (km/s) / (Mpc/h)
                  = dt_physical [s] * h_little / KM_PER_MPC

    The render_image.c code does:
        x += interp_frac * snap_dt * v
    so snap_dt needs to absorb the 1/a factor too... but a changes across
    the interval.  We use the mean a = (a1+a2)/2 as a reasonable approximation
    for the interpolation midpoint.  For interp_frac < 0.5 the error is small.

    Returns
    -------
    snap_dt : float
        In units of (Mpc/h) / (km/s).  Pass this directly to render_image.
    dt_gyr : float
        Physical time interval in Gyr (for sanity checking).
    """
    dt_s = comoving_time_integral(a1, a2, H0, Omega_m, Omega_L, n_steps)

    # Convert physical time to snap_dt in (Mpc/h)/(km/s)
    # snap_dt [Mpc/h / km/s] = dt [s] * h_little / KM_PER_MPC
    snap_dt = dt_s * h_little / KM_PER_MPC

    # For reference: Gyr
    SEC_PER_GYR = 3.15576e16
    dt_gyr = dt_s / SEC_PER_GYR

    return snap_dt, dt_gyr


def read_hdf5_params(filename, convention='GADGET4'):
    """Read a, H0, Omega_m, Omega_L, h from an HDF5 snapshot header."""
    try:
        import h5py
    except ImportError:
        sys.exit("h5py not installed. Run: pip install h5py")
    
    with h5py.File(filename, "r") as f:
        if convention == 'SWIFT':
            a = f['Header'].attrs['Scale-factor'] 
            H0 = f['Cosmology'].attrs['H0 [internal units]']
            h = f['Cosmology'].attrs['h'] 
            omega_bar   = f['Cosmology'].attrs['Omega_b']
            omega_dm    = f['Cosmology'].attrs['Omega_cdm']
            Omega_m    = omega_bar + omega_dm
            Omega_L    = f['Cosmology'].attrs['Omega_lambda'] 
        elif convention in ['GADGET4', 'AREPO']:
            a = f['Header'].attrs['Time'] 
            H0 = f['Parameters'].attrs['HubbleParam'] * 100.0
            h = f['Parameters'].attrs['HubbleParam'] 
            Omega_m = f['Parameters'].attrs['Omega0'] 
            Omega_L = f['Parameters'].attrs['OmegaLambda']
        else:
            a = f['Header'].attrs['Time'] 
            H0 = f['Header'].attrs['HubbleParam'] * 100.0
            h = f['Header'].attrs['HubbleParam'] 
            Omega_m = f['Header'].attrs['Omega0'] 
            Omega_L = f['Header'].attrs['OmegaLambda']
    return a, H0, Omega_m, Omega_L, h

def main():
    parser = argparse.ArgumentParser(
        description="Compute snap_dt for render_image velocity interpolation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    # HDF5 auto-read
    parser.add_argument("--hdf5",      metavar="FILE",
                        help="First snapshot (reads a, cosmology from header)")
    parser.add_argument("--hdf5_next", metavar="FILE",
                        help="Second snapshot (reads a from header)")
    parser.add_argument("--hdf5_convention", choices=['GADGET4', 'AREPO', 'SWIFT'], default='GADGET4',
                        help="HDF5 convention for reading parameters (default: GADGET4)")
    
    # Manual cosmology
    parser.add_argument("--snap_a1",  type=float, help="Scale factor of snapshot 1")
    parser.add_argument("--snap_a2",  type=float, help="Scale factor of snapshot 2")
    parser.add_argument("--snap_z1",  type=float, help="Redshift of snapshot 1 (alternative to --snap_a1)")
    parser.add_argument("--snap_z2",  type=float, help="Redshift of snapshot 2 (alternative to --snap_a2)")
    parser.add_argument("--H0",       type=float, default=67.74,
                        help="Hubble constant in km/s/Mpc (default: 67.74)")
    parser.add_argument("--Omega_m",  type=float, default=0.309,
                        help="Matter density parameter (default: 0.309)")
    parser.add_argument("--Omega_L",  type=float, default=0.691,
                        help="Dark energy density parameter (default: 0.691)")
    parser.add_argument("--h",        type=float, default=None,
                        help="Little h = H0/100 (derived from --H0 if not given)")

    # Output detail
    parser.add_argument("--n_interp", type=int, default=10,
                        help="Number of interpolated sub-frames to show (default: 10)")

    args = parser.parse_args()

    # ── Resolve parameters ────────────────────────────────────────────────────
    if args.hdf5:
        print(args.hdf5_convention + " convention: reading parameters from HDF5 headers...")
        a1, H0, Omega_m, Omega_L, h = read_hdf5_params(args.hdf5, convention=args.hdf5_convention)
        print(f"Read from {args.hdf5}: a={a1.item():.5f}  H0={H0.item():.02f}  "
              f"Omega_m={Omega_m.item():.4f}  Omega_L={Omega_L.item():.4f}  h={h.item():.4f}")
        if args.hdf5_next:
            a2, *_ = read_hdf5_params(args.hdf5_next, convention=args.hdf5_convention)
            print(f"Read from {args.hdf5_next}: a={a2.item():.5f}")
        else:
            sys.exit("--hdf5_next required when --hdf5 is given")
    else:
        H0      = args.H0
        Omega_m = args.Omega_m
        Omega_L = args.Omega_L
        h       = args.h if args.h is not None else H0 / 100.0

        a1 = args.snap_a1 if args.snap_a1 else (1.0 / (1.0 + args.snap_z1)) if args.snap_z1 else None
        a2 = args.snap_a2 if args.snap_a2 else (1.0 / (1.0 + args.snap_z2)) if args.snap_z2 else None

        if a1 is None or a2 is None:
            parser.print_help()
            sys.exit("\nError: provide --snap_a1/a2, --snap_z1/z2, or --hdf5/hdf5_next")

    # ── Compute ───────────────────────────────────────────────────────────────
    snap_dt, dt_gyr = snap_dt_gadget(a1, a2, H0, Omega_m, Omega_L, h)

    a_mid = 0.5 * (a1 + a2)
    z1 = 1.0/a1 - 1.0
    z2 = 1.0/a2 - 1.0

    print()
    print("=" * 60)
    print(f"  Snapshot interval:  a = {a1.item():.5f} to {a2.item():.5f}")
    print(f"                      z = {z1.item():.4f} to {z2.item():.4f}")
    print(f"  Physical dt       = {dt_gyr.item():.4f} Gyr")
    print(f"  snap_dt           = {snap_dt.item():.6g}  (Mpc/h) / (km/s)")
    print()
    print("  Pass to render_image:")
    print(f"    -snap_dt {snap_dt.item():.6g}")
    print("=" * 60)
    print()

    # ── Sub-frame table ───────────────────────────────────────────────────────
    print(f"  Sub-frame interpolation ({args.n_interp} frames):")
    print(f"  {'Frame':>6}  {'interp_frac':>12}  {'Δx / v  [Mpc/h / (km/s)]':>26}")
    print("  " + "-" * 50)
    for k in range(args.n_interp + 1):
        frac = k / args.n_interp
        dx_per_v = frac * snap_dt
        print(f"  {k:>6}  {frac:>12.3f}  {dx_per_v:>26.6g}")
    print()

    # ── Sanity check: typical displacement ───────────────────────────────────
    # Typical dark matter peculiar velocity ~200 km/s
    v_typical = 200.0   # km/s
    dx_typical = v_typical * snap_dt   # Mpc/h over full interval
    dx_pix_768 = dx_typical / (1.0)    # assume lbox=1 Mpc/h → pixels at 768px
    print(f"  Sanity: v={v_typical:.0f} km/s moves {dx_typical:.4f} Mpc/h "
          f"over one snapshot interval")
    print(f"  At lbox=1 Mpc/h, 768px: that's {dx_typical*768:.1f} pixels per full step")
    print()

    # ── Shell command ──────────────────────────────────────────────────────────
    print("  Example render_image command (10 sub-frames, fixed levels):")
    print()
    print("    N=10")
    print("    for i in $(seq 0 $N); do")
    print(f"      FRAC=$(python3 -c \"print(${{i}} / $N)\")")
    print("      ./render_image -input snapshot_A \\")
    print("          -output frame_${i} -isHDF5 \\")
    print("          -xc <xc> -yc <yc> -zc <zc> -lbox <lbox> \\")
    print("          -dark_matter -colormap magma \\")
    print("          -no_auto_levels -vmin <vmin> -vmax <vmax> \\")
    print(f"          -interp_frac $FRAC -snap_dt {snap_dt:.6g}")
    print("    done")
    print()


if __name__ == "__main__":
    main()
