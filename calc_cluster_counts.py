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


Format of data file "dndm.dat"
------------------------------
z1  mass1  density1
z1  mass2  density2
z1  ..     density3
z2  mass1  density4
z2  mass2  density5
z2  ..     density6
...
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
    print("FoV coverage: %sx%s [deg^2] = %.f [sr]" %
          (args.width, args.height, coverage), file=sys.stderr)
    Mmin_halo = args.Mmin * args.fraction
    print("Halo minimum mass: %g * %s = %g [Msun]" %
          (args.Mmin, args.fraction, Mmin_halo), file=sys.stderr)

    dndmdata = np.loadtxt(args.data)
    print("Loaded dndM data from file: %s" % args.data, file=sys.stderr)

    redshifts = set(dndmdata[:, 0])
    redshifts = np.array(list(redshifts))
    redshifts.sort()
    dz = redshifts[1] - redshifts[0]  # even-spacing

    halomass = dndmdata[:, 1]
    midx = (halomass >= Mmin_halo)
    mass1 = halomass[midx]
    midx2 = np.concatenate([midx[1:], [False]])
    mass2 = halomass[midx2]
    dmass = mass1 - mass2

    redshifts = dndmdata[midx, 0]
    density = dndmdata[midx, 2]

    counts = 0
    for z, dM, den in zip(redshifts, dmass, density):
        dVc = cosmo.differential_comoving_volume(z).value  # [Mpc^3/sr]
        counts += (den * dVc*dz*coverage * dM)
    print("Total number of clusters (M > %g [Msun]):" % Mmin, file=sys.stderr)
    print(int(counts))


if __name__ == "__main__":
    main()
