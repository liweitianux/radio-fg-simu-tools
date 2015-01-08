#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Make image of Galactic *synchrotron* radiation at given frequency.
#
# Weitian LI <liweitianux@gmail.com>
# 2015/01/07
#
# ChangeLogs:
# 2015/01/08, Weitian LI
#   * Added more detail 'USAGE', and improved usage()
#   * added argument 'verbose' for mkimg_galsyn()
#

from astropy.io import fits
import numpy as np

import sys
import getopt
import datetime


# Frequency of the raw input image (MHz)
raw_freq = 408


USAGE = """Usage:
    %(prog)s [ -h -C -v ] -f freq_MHz -i indexfile -r rawfile -o outfile

Required arguments:
    -f, --freq
        frequency (MHz) at which the radiation image to be made
    -i, --index
        index FITS file containing the spectral index values for each pixel
    -r, --raw
        raw FITS file of the %(raw_freq)s MHz radiation image
    -o, --out
        output FITS file of the radiation image at given frequency

Optional arguments:
    -h, --help
        print this usage
    -C, --clobber
        overwrite output file if already exists
    -v, --verbose
        show verbose information
""" % {'prog': sys.argv[0], 'raw_freq': raw_freq}


def mkimg_galsyn(freq, indexfile, rawfile, outfile,
        clobber=False, verbose=False):
    """
    Make image of Galactic synchrotron at given frequency.
    """
    if verbose:
        print("Openning %s and %s ..." % (indexfile, rawfile))
    index_fits  = fits.open(indexfile)
    raw_fits    = fits.open(rawfile)
    index_shape = index_fits[0].shape
    raw_shape   = raw_fits[0].shape
    index_data  = index_fits[0].data
    raw_data    = raw_fits[0].data
    assert (index_shape == raw_shape), \
            "Different shape of '%s' and '%s'!" % (indexfile, rawfile)
    # outfile shape
    rows = index_shape[0]
    cols = index_shape[1]
    out_data = np.zeros((rows, cols))
    # calculate outfile data
    if verbose:
        print("Calculating output image data ...")
    for i in range(rows):
        for j in range(cols):
            index = index_data[i, j]
            out_data[i, j] = raw_data[i, j] * np.power(freq/raw_freq, index)
    # copy header of rawfile to outfile
    out_header = raw_fits[0].header.copy(strip=True)
    out_header['DATE'] = (
            datetime.datetime.today().strftime('%Y-%m-%dT%H:%M:%S'),
            "Created with astropy.io.fits.")
    # write data and header to fits file
    if verbose:
        print("Writing data to %s ..." % outfile)
    out_fits = fits.PrimaryHDU(data=out_data, header=out_header)
    out_fits.writeto(outfile, clobber=clobber, checksum=True)
    # close fits file
    index_fits.close()
    raw_fits.close()


def usage():
    print(USAGE)


def main():
    """
    Process command line arguments, and make image.
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "Chf:i:r:o:v",
                ["clobber", "help", "freq=", "index=",
                 "raw=", "out=", "verbose"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    verbose = False
    clobber = False
    for opt, arg in opts:
        if opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-C", "--clobber"):
            clobber = True
        elif opt in ("-f", "--freq"):
            freq = float(arg)
        elif opt in ("-i", "--index"):
            indexfile = arg
        elif opt in ("-r", "--raw"):
            rawfile = arg
        elif opt in ("-o", "--out"):
            outfile = arg
        else:
            assert False, "unhandled option"

    if verbose:
        print("freq      = %s MHz" % freq)
        print("indexfile = %s" % indexfile)
        print("rawfile   = %s" % rawfile)
        print("outfile   = %s" % outfile)
        print("clobber   = %s" % clobber)

    # Make image
    mkimg_galsyn(freq=freq, indexfile=indexfile, rawfile=rawfile,
            outfile=outfile, clobber=clobber, verbose=verbose)


if __name__ == "__main__":
    main()

