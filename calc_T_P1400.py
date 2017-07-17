#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Based on 'calc_T_P1400.cc'
#
# Aaron LI
# 2015/04/16
#

"""
This tool calculate the X-ray temperature and power density at 1.4GHz
for each cluster.
"""


from astropy.cosmology import FlatLambdaCDM

import sys
import os
import re
import math
import datetime


## constants
H0 = 71 # Hubble constant at z=0 (km/s/Mpc)
Om0 = 0.27 # Omega matter at z=0
Tcmb0 = 2.725 # Temperature of the CMB at z=0 (default: 2.725 K)

## unit conversion
erg2Joule = 1.0e-7 # 1 erg = 1.0e-7 Joule

## cosmology calculator
cosmo = FlatLambdaCDM( H0=H0, Om0=Om0, Tcmb0=Tcmb0 )


def main():
    if len( sys.argv ) != 3:
        print( "Usage:\n    %s <in: z_m.dat> <out: z_m_T_P1400.dat>" %
                os.path.basename(sys.argv[0]), file=sys.stderr )
        sys.exit( 1 )

    infile = open( sys.argv[1], "r" )
    outfile = open( sys.argv[2], "w" )
    # Write file header
    outfile.write( "# Created: %(prog)s, %(time)s\n" %
            { 'prog': os.path.basename(sys.argv[0]),
              'time': datetime.datetime.now().isoformat() } )
    outfile.write( "# z   m(Msun)   Tx(keV)   P1400(W/Hz)\n" )

    for line in infile:
        if re.match( r"^\s*#", line ):
            continue
        z, m = map( float, line.split() )
        # Calculate X-ray temperature of the cluster with the
        # mass-temperature relation: (Arnaud et al. 2005)
        # m: unit Msun; Tx: unit keV
        Tx = 5.0 * math.pow( (m * cosmo.efunc(z) / 5.34e14), 1.0/1.72 )
        # The parameter values of this P_{1.4GHz} - Tx relation
        # was fitted by self with the data from Cassano et al. (2006):
        # Tx: unit keV; P1400: unit erg/s/Hz -> W/Hz
        # Formula: equation (5) of Wang et al. (2010)
        #P1400 = math.pow(T,2.88476)*1.072655e22*70*70/71./71.
        P1400 = 1.04e29*erg2Joule * math.pow(Tx, 2.88)
        outfile.write( "%.10lg %.10lg %.10lg %s\n" % (z, m, Tx, P1400) )

    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()

