#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Make image of Galactic *synchrotron* radiation at given frequency.
#
# TODO:
# * expand script description
# * describe the theory/calculation
# * provide references
#
# Aaron LI <aaronly.me@gmail.com>
# 2015/01/07
#
# ChangeLogs:
# 2017-05-05, Aaron LI
#   * Add option "-s / --single" to save image in single-precision float
# 2015/03/30, Aaron LI
#   * Add header keywords "FREQ" & "COMP"
# 2015/03/28, Aaron LI
#   * update author name
#   * replace for loops with numpy array manipulation in the output image
#     calculation, which greatly improve the calculate speed.
#   * update parameter name ['index', 'raw', 'out'] to
#     ['indexfile', 'rawfile', 'outfile'].
# 2015/03/26, Aaron LI
#   * Rename 'mkimg_galsync()' to 'mkfits_galsync()', and split the
#     write file part into 'main()'
#   * Add 'cmd_list' to records the tool and its parameters
#   * Add 'mkcards_hist()' to generate history cards
#   * Add HISTORY of process tool and parameters to output FITS header
# 2015/01/08, Aaron LI
#   * Added more detail 'USAGE', and improved usage()
#   * added argument 'verbose' for mkimg_galsyn()
#

from astropy.io import fits
import numpy as np

import os
import sys
import getopt
import datetime


# Frequency of the raw input image (MHz)
raw_freq = 408


USAGE = """Usage:
    %(prog)s [ -h -C -s -v ] -f freq_MHz -i indexfile -r rawfile -o outfile

Required arguments:
    -f, --freq
        frequency (MHz) at which the radiation image to be made
    -i, --indexfile
        index FITS file containing the spectral index values for each pixel
    -r, --rawfile
        raw FITS file of the %(raw_freq)s MHz radiation image
    -o, --outfile
        output FITS file of the radiation image at given frequency

Optional arguments:
    -h, --help
        print this usage
    -C, --clobber
        overwrite output file if already exists
    -s, --single
        save image in single-precision float (instead of double-precision)
    -v, --verbose
        show verbose information
""" % {'prog': os.path.basename(sys.argv[0]), 'raw_freq': raw_freq}


def mkfits_galsyn(freq, indexfile, rawfile, float32=False, verbose=False):
    """
    Return FITS of Galactic synchrotron image at given frequency,
    with header copyed from the rawfile.
    """
    if verbose:
        print("Openning %s and %s ..." % (indexfile, rawfile))
    index_fits  = fits.open(indexfile)
    raw_fits    = fits.open(rawfile)
    index_data  = index_fits[0].data.astype(np.float64)
    raw_data    = raw_fits[0].data.astype(np.float64)
    index_shape = index_fits[0].shape
    raw_shape   = raw_fits[0].shape
    assert (index_shape == raw_shape), \
            "Different shape of '%s' and '%s'!" % (indexfile, rawfile)
    # calculate outfile data
    if verbose:
        print("Calculating output image data ...")
    out_data = raw_data * np.power(freq/raw_freq, index_data)
    if float32:
        out_data = out_data.astype(np.float32)
    # copy header of rawfile to outfile
    out_header = raw_fits[0].header.copy(strip=True)
    # Add meta information of the new fits to header
    out_header.set("FREQ", freq, "MHz")
    out_header.set("COMP", "Galactic_synchrotron", "Radiation component")
    # close fits file
    index_fits.close()
    raw_fits.close()
    # create the FITS object including the image data and header
    out_fits = fits.PrimaryHDU(data=out_data, header=out_header)
    return out_fits


def mkcards_hist(cmd_list):
    """
    Return a list of HISTORY Card, which record the tool and parameters
    used in this process.
    """
    tool_card = fits.Card("HISTORY",
            "TOOL: %(tool)s (%(time)s)" % {
                'tool': cmd_list[0],
                'time': datetime.datetime.today().strftime("%Y-%m-%dT%H:%M:%S") },
            "by Weitian LI (c) 2015")
    parm = "PARM: " + " ".join(cmd_list[1:])
    parm_card = fits.Card("HISTORY", parm)
    return [tool_card, parm_card]


def usage():
    print(USAGE)


def main():
    """
    Process command line arguments, and make image.
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "Chf:i:r:o:sv",
                ["clobber", "help", "freq=", "indexfile=",
                 "rawfile=", "outfile=", "single", "verbose"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    verbose = False
    clobber = False
    # Records this tool and its options/arguments,
    # and used to write history into FITS header.
    cmd_list = [ os.path.basename(sys.argv[0]) ]
    for opt, arg in opts:
        if opt in ("-v", "--verbose"):
            verbose = True
            cmd_list.append("--verbose")
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-C", "--clobber"):
            clobber = True
            cmd_list.append("--clobber")
        elif opt in ("-f", "--freq"):
            freq = float(arg)
            cmd_list.append("--freq=%s" % freq)
        elif opt in ("-i", "--indexfile"):
            indexfile = arg
            cmd_list.append("--indexfile=%s" % indexfile)
        elif opt in ("-r", "--rawfile"):
            rawfile = arg
            cmd_list.append("--rawfile=%s" % rawfile)
        elif opt in ("-o", "--outfile"):
            outfile = arg
            cmd_list.append("--outfile=%s" % outfile)
        elif opt in ("-s", "--single"):
            float32 = True
            cmd_list.append("--single")
        else:
            assert False, "unhandled option"

    if verbose:
        print("freq      = %s MHz" % freq)
        print("indexfile = %s" % indexfile)
        print("rawfile   = %s" % rawfile)
        print("outfile   = %s" % outfile)
        print("clobber   = %s" % clobber)
        print("float32   = %s" % float32)

    # Create the FITS object for the radiation image
    out_fits = mkfits_galsyn(freq=freq, indexfile=indexfile,
                             rawfile=rawfile, float32=float32,
                             verbose=verbose)

    # Update fits header to record this tool and parameters.
    hist_cards = mkcards_hist(cmd_list)
    out_fits.header.extend(hist_cards)

    # Write FITS object into output file.
    if verbose:
        print("Writing data to %s ..." % outfile)
    out_fits.writeto(outfile, clobber=clobber, checksum=True)


if __name__ == "__main__":
    main()

