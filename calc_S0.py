#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Based on 'calc_s0.cc' & 'calc_F1400_rc.py'
#
# Aaron LI
# 2015/04/21
#

"""
This tool calculate the S0 parameter of the beta-model.
"""


from astropy.cosmology import FlatLambdaCDM
from scipy import integrate

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
kg = 1000 # CGS unit system
M_sun = 1.989e30 * kg # mass of sun (unit: g)

## unit conversion
erg2Joule = 1.0e-7 # 1 erg = 1.0e-7 Joule

## cosmology calculator
cosmo = FlatLambdaCDM( H0=H0, Om0=Om0, Tcmb0=Tcmb0 )


def rad2arcmin( rad ):
    """
    Convert value in unit 'rad' to 'arcmin'.
    """
    return rad / math.pi * 180.0 * 60.0


def f_beta( r, *args ):
    """
    Integrand function to calculate the S0 parameter of beta-model.
    """
    rc = args[0]
    beta = args[1]
    return 2.0*math.pi * r * math.pow( 1.0 + (r/rc)**2, 0.5 - 3.0*beta )


def main():
    if len( sys.argv ) != 3:
        print( "Usage:\n    %s " % os.path.basename(sys.argv[0]) +
                "<in: z_m_F${nu}_rc_T.dat> <out: F${nu}_S0_r200_rc_beta.dat>",
                file=sys.stderr )
        sys.exit( 1 )

    infile = open( sys.argv[1], "r" )
    outfile = open( sys.argv[2], "w" )
    # Write file header
    outfile.write( "# Created: %(prog)s, %(time)s\n" %
            { 'prog': os.path.basename(sys.argv[0]),
              'time': datetime.datetime.now().isoformat() } )
    outfile.write( "# F${nu}(erg/s/cm^2)   S0(???)   r200(arcmin)   " \
            + "rc(arcmin)   beta\n" )

    for line in infile:
        if re.match( r"^\s*#", line ):
            continue
        z, m, flux, rc_arcmin, Tx = map( float, line.split() )
        #
        rho_c = cosmo.critical_density( z ).value # unit: g/cm^3
        delta_c = 200 # calculate r200
        # Calculate the viral radius R200
        # (4*pi/3) * r200^3 * (200*rho_c) = m * M_sun
        r200_cm = math.pow( 3.0*m*M_sun / (4.0*math.pi * delta_c*rho_c), 1.0/3.0 )
        da_cm = cosmo.angular_diameter_distance( z ).value * Mpc
        # convert r200_cm from unit cm to arcmin & Mpc
        r200_arcmin = rad2arcmin( r200_cm / da_cm )
        r200_Mpc = r200_cm / Mpc
        # Calculate the beta parameter of beta-model
        # Reference: equation (10) of Wang et al. (2010)
        beta = 8.85e-15 * m * (1.0 + (r200_arcmin/rc_arcmin)**2) \
                / ((r200_arcmin/rc_arcmin)**2 * Tx * r200_Mpc)
        value = integrate.quad( f_beta, 0.0, r200_arcmin,
                args=(rc_arcmin, beta) )
        S0 = flux / value[0]
        #
        outfile.write( "%.10lg %.10lg %.10lg %.10lg %.10lg\n"
                % (flux, S0, r200_arcmin, rc_arcmin, beta) )

    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()

