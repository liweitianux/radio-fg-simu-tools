#!/usr/bin/env python3
#
# Copyright (c) 2015,2017 Aaron LI
#
# Based on 'calc_totalnum.cc'
#
# 2015/04/16
# Updated: 2017-07-17
#

"""
This tool calculate the total number of clusters within the specified
sky coverage.


Format of data file "dndmdata.txt"
------------------------------
z1  mass1  density1
z1  mass2  density2
z1  ..     density3
z2  mass1  density4
z2  mass2  density5
z2  ..     density6
...

where,
mass has unit [Msun],
density has unit [number]/[Msun]/dVc, and dVc is differential comoving
volume with dimension [Mpc^3]/[unit redshift]/[sr]
"""


import sys
import argparse

import numpy as np
import astropy.units as au
from astropy.cosmology import FlatLambdaCDM


# Cosmology
H0 = 71.0  # Hubble constant at z=0 [km/s/Mpc]
Om0 = 0.27  # Critical matter density at z=0
cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)


def load_dndmdata(filename):
    """
    Load dndM data and reformat into a 2D density grid together with
    redshifts and masses vectors.

    Data Description
    ----------------
    * Redshifts: 0.0 -> 3.02, even-spacing, step 0.02
    * Halo mass: 1e12 -> 9.12e15 [Msun], log-even (dark matter)
    * Density: number/dVc/dM
      where,
      - dVc: differential comvoing volume, Mpc^3/sr/[unit redshift]
    """
    data = np.loadtxt(filename)
    redshifts = data[:, 0]
    masses = data[:, 1]
    density = data[:, 2]

    redshifts = np.array(list(set(redshifts)))
    redshifts.sort()
    masses = np.array(list(set(masses)))
    masses.sort()
    density = density.reshape((len(redshifts), len(masses)))

    return (redshifts, masses, density)


def delta(x, logeven=False):
    """
    Calculate the delta values for each element of a vector,
    assuming they are evenly or log-evenly distributed,
    with extrapolating.
    """
    x = np.asarray(x)
    if logeven:
        x = np.log(x)
    step = x[1] - x[0]
    x1 = np.concatenate([[x[0]-step], x[:-1]])
    x2 = np.concatenate([x[1:], [x[-1]+step]])
    dx = (x2 - x1) * 0.5
    if logeven:
        dx = np.exp(dx)
    return dx


def calc_number_grid(density, redshifts, masses):
    """
    Calculate the number distribution w.r.t. redshift, mass, and unit
    coverage [sr] from the density distribution.
    """
    dz = delta(redshifts)
    dM = delta(masses)
    dMgrip, dzgrip = np.meshgrid(dM, dz)
    Mgrip, zgrip = np.meshgrid(masses, redshifts)
    dVcgrip = cosmo.differential_comoving_volume(zgrip).value  # [Mpc^3/sr]
    numgrid = density * dVcgrip * dzgrip * dMgrip
    return numgrid


def calc_cluster_counts(numgrid, masses, Mmin, coverage):
    """
    Calculate the total number of clusters (>= minimum mass) within the
    FoV coverage according to the number density distribution (e.g.,
    predicted by the Press-Schechter mass function)
    """
    midx = masses >= Mmin
    counts = np.sum(numgrid[:, midx]) * coverage
    return counts


def main():
    # Default parameters
    fov_width = 10.0  # Width of FoV [deg]
    fov_height = fov_width  # Height of FoV [deg]
    dm_fraction = 0.83  # dark matter fraction in clusters
    Mmin = 2e14  # minimum (total) mass of clusters

    parser = argparse.ArgumentParser(
        description="Calculate total cluster number (>Mmin) within a FoV")
    parser.add_argument("-W", "--width", dest="width", type=float,
                        default=fov_width,
                        help="FoV width [deg] (default: %s)" % fov_width)
    parser.add_argument("-H", "--height", dest="height", type=float,
                        default=fov_height,
                        help="FoV height [deg] (default: %s)" % fov_height)
    parser.add_argument("-f", "--fraction", dest="fraction", type=float,
                        default=dm_fraction,
                        help="dark matter fraction in clusters " +
                        "(default: %s)" % dm_fraction)
    parser.add_argument("-m", "--min-mass", dest="Mmin", type=float,
                        default=Mmin,
                        help="Minimum (total) mass of the clusters " +
                        "(default: %s [Msun])" % Mmin)
    parser.add_argument("-d", "--data", dest="data", required=True,
                        help="dndM data file generated from Press-Schechter " +
                        "formalism (e.g., dndMdata.txt)")
    args = parser.parse_args()

    coverage_deg2 = args.width * args.height  # [deg^2]
    coverage = coverage_deg2 * au.deg.to(au.rad)**2  # [rad^2] = [sr]
    print("FoV coverage: %sx%s [deg^2] = %g [sr]" %
          (args.width, args.height, coverage), file=sys.stderr)
    Mmin_halo = args.Mmin * args.fraction
    print("Halo minimum mass: %g * %s = %g [Msun]" %
          (args.Mmin, args.fraction, Mmin_halo), file=sys.stderr)

    redshifts, masses, density = load_dndmdata(args.data)
    numgrid = calc_number_grid(density, redshifts=redshifts, masses=masses)
    counts = calc_cluster_counts(numgrid, masses=masses, Mmin=Mmin_halo,
                                 coverage=coverage)

    print("Total number of clusters (M > %g [Msun]):" % Mmin, file=sys.stderr)
    print(int(counts))


if __name__ == "__main__":
    main()
