#!/usr/bin/env python3
#
# Copyright (c) 2015,2017 Aaron LI
#
# Based on `calc_totalnum.cc` & `z_m_simulation.cc`
#
# 2015/04/16
# Updated: 2017-07-19
#

"""
Calculate the total number of clusters within the specified FoV coverage
according to the supplied dark matter density distribution obtained by
Press-Schechter formalism, then sampling the density distribution to
generate random (z, M) pairs for each cluster using the "rejection algorithm".


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


import os
import sys
import argparse
import random

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
    return int(np.round(counts))


def bilinear_interpolation(x, y, p11, p12, p21, p22):
    """
    Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

    p11 --+-- p12
     |    |    |
     +--(x,y)--+
     |    |    |
    p21 --+-- p22

    Credit
    ------
    [1] http://en.wikipedia.org/wiki/Bilinear_interpolation
    [2] https://stackoverflow.com/a/8662355/4856091
    """
    # Sort points by X then Y
    points = sorted([p11, p12, p21, p22])
    x1, y1, q11 = points[0]
    _x1, y2, q12 = points[1]
    x2, _y1, q21 = points[2]
    _x2, _y2, q22 = points[3]

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError("points do not form a rectangle")
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError("(x, y) not within the rectangle")
    a1 = (q11 * (x2 - x) * (y2 - y) + q21 * (x - x1) * (y2 - y) +
          q12 * (x2 - x) * (y - y1) + q22 * (x - x1) * (y - y1))
    a2 = (x2 - x1) * (y2 - y1)
    q = a1 / a2
    return q


def sample_z_m(num, numgrid, redshifts, masses, Mmin):
    """
    Randomly generate the requested number of pairs of (z, M) following
    the specified number distribution.
    """
    zmin = redshifts.min()
    zmax = redshifts.max()
    Mmax = masses.max()
    midx = masses >= Mmin
    numgrid2 = numgrid[:, midx]
    NM = numgrid2.max()
    z_list = []
    M_list = []
    i = 0
    while i < num:
        z = random.uniform(zmin, zmax)
        M = random.uniform(Mmin, Mmax)
        r = random.random()
        zi1 = (redshifts < z).sum()
        zi2 = zi1 - 1
        if zi2 < 0:
            zi2 += 1
            zi1 += 1
        Mi1 = (masses < M).sum()
        Mi2 = Mi1 - 1
        if Mi2 < 0:
            Mi2 += 1
            Mi1 += 1
        N = bilinear_interpolation(
            z, np.log(M),
            p11=(redshifts[zi1], np.log(masses[Mi1]), numgrid[zi1, Mi1]),
            p12=(redshifts[zi1], np.log(masses[Mi2]), numgrid[zi1, Mi2]),
            p21=(redshifts[zi2], np.log(masses[Mi1]), numgrid[zi2, Mi1]),
            p22=(redshifts[zi2], np.log(masses[Mi2]), numgrid[zi2, Mi2]))
        if r < N/NM:
            print(i+1, "/", num, "...")
            z_list.append(z)
            M_list.append(M)
            i += 1
    return (np.array(z_list), np.array(M_list))


def main():
    # Default parameters
    fov_width = 10.0  # Width of FoV [deg]
    fov_height = fov_width  # Height of FoV [deg]
    dm_fraction = 0.80  # dark matter fraction in clusters
    Mmin = 2e14  # minimum (total) mass of clusters

    parser = argparse.ArgumentParser(
        description="Calculate total cluster number (>Mmin) within a FoV")
    parser.add_argument("-C", "--clobber", dest="clobber",
                        action="store_true",
                        help="overwrite existing file")
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
                        "(default: %g [Msun])" % Mmin)
    parser.add_argument("-d", "--data", dest="data", required=True,
                        help="dndM data file generated from Press-Schechter " +
                        "formalism (e.g., dndMdata.txt)")
    parser.add_argument("-S", "--sample-z-m", dest="sample_zm",
                        action="store_true",
                        help="whether to sample the (z, M) pairs")
    parser.add_argument("--csv", dest="csv", action="store_true",
                        help="whether to store the sample data in CSV format")
    parser.add_argument("-o", "--outfile", dest="outfile",
                        help="output file to save the sampled (z,M) pairs")
    args = parser.parse_args()

    if args.outfile and os.path.exists(args.outfile) and not args.clobber:
        raise OSError("output file already exists: %s" % args.outfile)

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
    print(counts)

    if args.sample_zm:
        if args.outfile is None:
            raise ValueError("--outfile required!")
        print("Sampling (z, M) pairs ...")
        randz, randM = sample_z_m(num=counts, numgrid=numgrid,
                                  redshifts=redshifts,
                                  masses=masses, Mmin=Mmin_halo)
        randM /= args.fraction  # -> cluster mass [Msun]
        catalog = np.column_stack([randz, randM])
        header = ["dndM data: %s\n" % os.path.abspath(args.data),
                  "coverage: %s x %s [deg]\n" % (args.width, args.height),
                  "dark matter fraction: %s\n" % args.fraction,
                  "minimum cluster mass: %g [Msun]\n" % args.Mmin,
                  "cluster numbers: %d\n" % counts,
                  "----------------------------------------\n",
                  "redshift               ClusterMass[Msun]"]
        if args.csv:
            print("Save sampled (z, M) in CSV format ...", file=sys.stderr)
            import pandas as pd

            df = pd.DataFrame(catalog, columns=["z", "mass"])
            with open(args.outfile, "w") as fp:
                fp.write("# " + "# ".join(header) + "\n")
                df.to_csv(fp, index=False)
        else:
            np.savetxt(args.outfile, np.column_stack([randz, randM]),
                       header="".join(header))
        print("Sampled (z, M) data written to: %s" % args.outfile,
              file=sys.stderr)


if __name__ == "__main__":
    main()
