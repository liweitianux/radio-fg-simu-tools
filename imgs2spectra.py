#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Convert a list of FITS images of different frequencies to
# a list of spectra of different pixel position, and save as
# HDF5 format with 'h5py' package.
#
# Aaron LI
# 2015/03/31
#
# ChangeLogs:
# 2015/03/31, Aaron
#   * open FITS file when used, and close once used.
#

from astropy.io import fits
import numpy as np
import h5py

import os
import sys
import getopt
import datetime


USAGE = """Usage:
    %(prog)s [ -h -C -v -z ] -c component -l FITS_list -o outfile

Required arguments:
    -c, --component
        Component type of radiation to be processed
    -l, --list
        Text file containing the list of FITS files
        These file should be sorted according to there frequencies.
    -o, --outfile
        Output HDF5 file to save the spectra information
    -z, --gzip
        Enable "gzip" compression for HDF5 dataset

Optional arguments:
    -h, --help
        print this usage
    -C, --clobber
        overwrite output file if already exists
    -v, --verbose
        show verbose information
""" % { 'prog': os.path.basename(sys.argv[0]) }


def usage():
    print(USAGE)


def main():
    """
    Process command line arguments, and calculate output image,
    update header history, and write to file.
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], "c:Chl:o:vz",
                ["component=", "clobber", "help", "list=",
                 "outfile=", "verbose", "gzip"])
    except getopt.GetoptError as err:
        print(err)
        usage()
        sys.exit(2)
    verbose = False
    clobber = False
    compression = None
    # Records this tool and its options/arguments,
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
        elif opt in ("-c", "--component"):
            component = arg
            cmd_list.append("--component=%s" % component)
        elif opt in ("-l", "--list"):
            listfile = arg
            cmd_list.append("--list=%s" % listfile)
        elif opt in ("-o", "--outfile"):
            outfile = arg
            cmd_list.append("--outfile=%s" % outfile)
        elif opt in ("-z", "--gzip"):
            compression = "gzip"
            cmd_list.append("--gzip")
        else:
            assert False, "unhandled option"

    if verbose:
        print("component   = %s" % component)
        print("list        = %s" % listfile)
        print("outfile     = %s" % outfile)
        print("clobber     = %s" % clobber)
        print("compression = %s" % compression)

    # open listfile and strip the newline character
    fn_list = [line.strip() for line in open(listfile)]

    # get basic information
    freq_list = []
    freq_num = len(fn_list)
    fits_tmp = fits.open(fn_list[0])
    rows, columns = fits_tmp[0].data.shape

    # open hdf5 file to save data
    if verbose:
        print("Open HDF5 file ...")
    if clobber:
        try:
            os.remove(outfile)
        except OSError:
            pass
    h5outfile = h5py.File(outfile, "w")

    # create group for this component
    comp_grp = h5outfile.create_group(component.upper())
    # add meta information to attributes
    comp_grp.attrs.create("component", component.upper(),
            dtype=h5py.special_dtype(vlen=str))
    comp_grp.attrs.create("freq_num", freq_num)
    comp_grp.attrs.create("rows", rows)
    comp_grp.attrs.create("columns", columns)
    # add process history to attributes
    history = mk_history(cmd_list)
    if verbose:
        print(history)
    for item in history:
        comp_grp.attrs.create(item["key"], item["value"],
                dtype=h5py.special_dtype(vlen=str))

    # create dataset to record spectra for each pixel
    if verbose:
        print("Create the spectra dataset ...")
    spectra_ds = comp_grp.create_dataset("spectra",
            (rows, columns, freq_num), dtype="float64",
            compression=compression)

    # convert images of different frequncies into spectra of different pixels
    if verbose:
        print("Convert images into spectra of different pixels ...")
    for k in range(freq_num):
        # open fits file
        fitsobj = fits.open(fn_list[k])
        freq = fitsobj[0].header["FREQ"]
        freq_list.append(freq)
        if verbose:
            print("frequency: %s" % freq)
        # assign FITS data
        spectra_ds[:, :, k] = fitsobj[0].data
        # close this fits file
        fitsobj.close()

    # create dataset to record frequncy information
    if verbose:
        print("Create the frequencies dataset ...")
    freq_ds = comp_grp.create_dataset("frequencies",
            (freq_num,), dtype="float64")
    freq_ds[:] = np.array(freq_list)

    # close hdf5 file
    h5outfile.close()


def mk_history(cmd_list):
    """
    Return a list of dictionaries which record the process history,
    i.e., this tool and its parameters.
    """
    history = [ {
        "key": "TOOL",
        "value": "%(tool)s (%(time)s)" % {
                'tool': cmd_list[0],
                'time': datetime.datetime.today().isoformat() } } ]
    history.append({
        "key": "PARAMETERS",
        "value": " ".join(cmd_list[1:]) })
    return history


if __name__ == "__main__":
    main()

