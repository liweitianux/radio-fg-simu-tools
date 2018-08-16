#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Based on 'calc_F1400_rcDA.cc'
#
# Aaron LI
# 2015/04/21
#

"""
This tool calculate the 1.4GHz flux density and core radius parameter of
the beta-model.
"""


from astropy.cosmology import FlatLambdaCDM

import sys
import os
import re
import math
import datetime


## constants
H0 = 71 # Hubble constant at z=0 (km/s/Mpc)
h71 = 1 # We use H0=71
Om0 = 0.27 # Omega matter at z=0
Tcmb0 = 2.725 # Temperature of the CMB at z=0 (default: 2.725 K)
pc = 3.08567758e18 # 1 pc = ? cm
kpc = 1000 * pc # 1 kpc = ? cm
Mpc = 1000 * kpc # 1 Mpc = ? cm

## unit conversion
erg2Joule = 1.0e-7 # 1 erg = 1.0e-7 Joule

## cosmology calculator
cosmo = FlatLambdaCDM( H0=H0, Om0=Om0, Tcmb0=Tcmb0 )


def rad2arcmin( rad ):
    """
    Convert value in unit 'rad' to 'arcmin'.
    """
    return rad / math.pi * 180.0 * 60.0


def main():
    if len( sys.argv ) != 3:
        print( "Usage:\n    %s " % os.path.basename(sys.argv[0]) +
                "<in: z_P1400_Lrosat.dat> <out: z_F1400_rc.dat>",
                file=sys.stderr )
        sys.exit( 1 )

    infile = open( sys.argv[1], "r" )
    outfile = open( sys.argv[2], "w" )
    # Write file header
    outfile.write( "# Created: %(prog)s, %(time)s\n" %
            { 'prog': os.path.basename(sys.argv[0]),
              'time': datetime.datetime.now().isoformat() } )
    outfile.write( "# z   F1400(erg/s/cm^2)   rc(arcmin)\n" )

    for line in infile:
        if re.match( r"^\s*#", line ):
            continue
        z, m, Tx, P1400, Lrosat = map( float, line.split() )
        # Calculate the core radius parameter of the beta-model
        # L_rosat: unit erg/s; r_c: unit kpc
        # Reference: equation (8) of Wang et al. (2010)
        rc = 176.0 / h71 * math.pow( Lrosat/5.0e44, 0.2 ) # kpc
        da_cm = cosmo.angular_diameter_distance( z ).value * Mpc
        # convert rc from unit kpc to arcmin
        rc_arcmin = rad2arcmin( rc * kpc / da_cm )
        dl_cm = cosmo.luminosity_distance( z ).value * Mpc
        F1400 = P1400 * 1.0e7 / (4.0*math.pi * dl_cm*dl_cm)
        outfile.write( "%.10lg %.10lg %.10lg\n" % (z, F1400, rc_arcmin) )

    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()

