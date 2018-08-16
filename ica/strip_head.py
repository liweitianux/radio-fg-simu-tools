#!/usr/bin/env python
import pyfits
import sys

if len(sys.argv)!=3:
    print "Usage: strip_head.py input output"
    sys.exit(-1)

pyfits.HDUList(pyfits.PrimaryHDU(pyfits.open(sys.argv[1])[0].data)).writeto(sys.argv[2],clobber=True)
