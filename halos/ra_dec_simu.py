#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Based on 'sphxy.cc'
#
# Aaron LI
# 2015/04/17
#

"""
This tool simulates the (RA, DEC) coordinates for each cluster.
"""


import sys
import os
import math
import random


## constants
## simulation related
sky_width_deg  = 10  # simulated width of sky coverage (degree)
sky_height_deg = 10  # simulated height of sky coverage (degree)


def main():
    if len( sys.argv ) != 3:
        print( "Usage:\n    %s <number> <out: ra_dec.dat>" %
                os.path.basename(sys.argv[0]), file=sys.stderr )
        sys.exit( 1 )

    number = int( sys.argv[1] )
    outfile = open( sys.argv[2], "w" )

    for i in range( number ):
        x = random.uniform( -0.5*sky_width_deg, 0.5*sky_width_deg )
        y = random.uniform( -0.5*sky_height_deg, 0.5*sky_height_deg )
        ra = math.atan2(y, x)*180.0/math.pi + 180.0 # [-pi, pi] => [0, 360]
        dec = 90.0 - math.sqrt( x**2 + y**2 )
        outfile.write( "%f %f\n" % (ra, dec) )
    outfile.close()

if __name__ == "__main__":
    main()

