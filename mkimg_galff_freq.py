#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Make image of Galactic *free-free* radiation at given frequency.
#
# The re-combination line of H (especially Halpha line) traces the
# Galactic free-free radiation (Smoot 1998; Reynold & Haffner 2000).
# Finkbeiner (2003) provided the brightness map of Halpha (I_Halpha),
# and the Galactic free-free radiation at 30 GHz can be calculated as:
# T^{Gff}_{30GHz}(r) = 7.4e-6 * ( I_Halpha(r) / Rayleigh ) (K)
# where: 1 Rayleigh = 1e6 / (4*pi) photons/s/cm^2/sr
#
# The spectrum of Galactic free-free radiation can be generally described
# by a broken-powerlaw model, with a spectral index \alpha = 2.10 when
# frequency \nu <= 10 GHz (Shaver et al. 1999), and index \alpha = 2.15
# when \nu > 10 GHz (Bennett et al. 2003).
#      T^{Gff}_{\nu}(r) \proto \nu^{- \alpha}
#   => T^{Gff}_{\nu}(r) = T^{Gff}_{30GHz} * (\nu / 30GHz)^{- \alpha}
#
# XXX:
# * unit of input Halpha image: Rayleigh? or K?
#
# Aaron LI <aaronly.me@gmail.com>
# 2015/03/26
#
# ChangeLogs:
# 2015/03/30, Aaron LI
#   * Add header keywords "FREQ" & "COMP"
# 2015/03/28, Aaron LI
#   * update author name
#   * replace for loops with numpy array manipulation in the output image
#     calculation, which greatly improve the calculate speed.
#

from astropy.io import fits
import numpy as np

import os
import sys
import getopt
import datetime


# Spectral indexes of Galactic free-free radiation:
index_10ghz_above = 2.15 # If frequency is above 10 GHz
index_10ghz_below = 2.10 # If frequency is below 10 GHz


USAGE = """Usage:
    %(prog)s [ -h -C -v ] -f freq_MHz -i Halpha_img -o outfile

Required arguments:
    -f, --freq
        frequency (MHz) at which the radiation image to be made
    -i, --infile
        H_alpha image which used as the calculation template
    -o, --outfile
        output FITS file of the radiation image at given frequency

Optional arguments:
    -h, --help
        print this usage
    -C, --clobber
        overwrite output file if already exists
    -v, --verbose
        show verbose information
""" % { 'prog': os.path.basename(sys.argv[0]) }


def mkfits_galff(freq, infile, verbose=False):
    """
    Return FITS of Galactic free-free image at given frequency,
    with header copyed from the input file.
    """
    if verbose:
        print("Openning %s ..." % infile)
    halpha_fits  = fits.open(infile)
    # convert data type to "float64"
    halpha_data  = halpha_fits[0].data.astype(np.float64)
    # calculate outfile data
    if verbose:
        print("Calculating output image data ...")
    T_gff_30ghz = 7.4e-6 * halpha_data # unit: K
    if freq >= 10000:
        # freq >= 10 GHz: one single power law is sufficient
        out_data = T_gff_30ghz * np.power(freq/30000, -index_10ghz_above)
    else:
        # freq < 10 GHz: requires broken power law
        T_gff_10ghz = T_gff_30ghz * np.power(10000/30000, -index_10ghz_above)
        out_data = T_gff_10ghz * np.power(freq/10000, -index_10ghz_below)
    # Copy the header of infile to outfile
    out_header = halpha_fits[0].header.copy(strip=True)
    # Remove some unwanted keywords from header
    for key in ["CONTENT", "HDUNAME"]:
        if key in out_header:
            del out_header[key]
    # Add meta information of the new fits to header
    out_header.set("FREQ", freq, "MHz")
    out_header.set("COMP", "Galactic_free-free", "Radiation component")
    # close fits file
    halpha_fits.close()
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
    Process command line arguments, and calculate output image,
    update header history, and write to file.
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "Cf:hi:o:v",
                ["clobber", "freq=", "help", "infile=",
                 "outfile=", "verbose"])
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
            assert (freq > 0), "Specified frequency = %f <= 0!" % freq
            cmd_list.append("--freq=%s" % freq)
        elif opt in ("-i", "--infile"):
            infile = arg
            cmd_list.append("--infile=%s" % infile)
        elif opt in ("-o", "--outfile"):
            outfile = arg
            cmd_list.append("--outfile=%s" % outfile)
        else:
            assert False, "unhandled option"

    if verbose:
        print("freq    = %s MHz" % freq)
        print("infile  = %s" % infile)
        print("outfile = %s" % outfile)
        print("clobber = %s" % clobber)

    # Create the FITS object for the radiation image
    out_fits = mkfits_galff(freq=freq, infile=infile, verbose=verbose)

    # Update fits header to record this tool and parameters.
    hist_cards = mkcards_hist(cmd_list)
    out_fits.header.extend(hist_cards)

    # Write FITS object into output file.
    if verbose:
        print("Writing data to %s ..." % outfile)
    out_fits.writeto(outfile, clobber=clobber, checksum=True)


if __name__ == "__main__":
    main()

