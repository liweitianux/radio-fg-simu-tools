#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Based on 'calc_totalnum.cc'
#
# Aaron LI
# 2015/04/16
#

"""
This tool calculate the total number of clusters within the specified
sky coverage.
"""


from astropy.cosmology import FlatLambdaCDM

import sys
import os
import re
import math


## constants
H0 = 71 # Hubble constant at z=0 (km/s/Mpc)
Om0 = 0.27 # Omega matter at z=0
Tcmb0 = 2.725 # Temperature of the CMB at z=0 (default: 2.725 K)
cm = 1 # CGS unit system
s = 1 # CGS unit system
km = 1000 * 100 # -> cm
c = 299792458 * 100 # speed of light (cm/s)
## simulation related
sky_width_deg  = 10  # simulated width of sky coverage (degree)
sky_height_deg = 10  # simulated height of sky coverage (degree)

## cosmology calculator
cosmo = FlatLambdaCDM( H0=H0, Om0=Om0, Tcmb0=Tcmb0 )


def main():
    if len( sys.argv ) != 3:
        print( "Usage:\n    %s <mass_limit> <datafile>" %
                os.path.basename(sys.argv[0]), file=sys.stderr )
        sys.exit( 1 )

    mass_limit = float( sys.argv[1] )

    ## Read z, mass, flux data from file
    ## file format:
    ##     z1  m1  fx
    ##     z1  m2  fx
    ##     z1  ..  fx
    ##     z2  m1  fx
    ##     z2  m2  fx
    ##     z2  ..  fx
    ##     ..
    z_list = []
    m_list = []
    f_list = []
    print( "Reading z, mass, flux data from file ...", file=sys.stderr )
    with open( sys.argv[2], "r" ) as datafile:
        for line in datafile:
            if re.match( r"^\s*#", line ):
                continue
            z, m, f = map( float, line.split() )
            if z not in z_list:
                z_list.append( z )
            if m not in m_list:
                m_list.append( m )
            f_list.append( f )
    print( "Data loaded.", file=sys.stderr )
    print( "z_num: %d; m_num: %d; f_num: %d" %
            (len(z_list), len(m_list), len(f_list)), file=sys.stderr )

    number = 0  # cluster numbers
    sky_coverage = (sky_width_deg * math.pi / 180) * \
            (sky_height_deg * math.pi / 180)
    for i in range( len(z_list)-1 ):
        for j in range( 1, len(m_list) ):
            if m_list[j] >= mass_limit:
                z = z_list[i]
                m = m_list[j]
                try:
                    f = f_list[ len(m_list)*i + j ]
                except IndexError:
                    print( "IndexError: index=%d out of range!" %
                            (len(z_list)*i + j), file=sys.stderr )
                    sys.exit(10)

                delta_z = z_list[i+1] - z
                delta_m = m - m_list[j-1]
                da_Mpc = cosmo.angular_diameter_distance( z ).value
                Hz = cosmo.H( z ).value
                number += delta_z * delta_m * f * (c/km) \
                        / Hz * da_Mpc**2 * (1+z)**2 * sky_coverage
    print( "cluster_number: %d" % int(number) )


if __name__ == "__main__":
    main()

