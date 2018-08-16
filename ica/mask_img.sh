#!/usr/bin/bash

dmimgcalc infile=Timg_bin_67_sm.fits infile2=none outfile=mask.fits op="imgout=img1-img1+1" clobber=yes


for i in 61 63 64 65 66 67 69
do
    ../ica_tools/mask_points.py Timg_bin_${i}_sm.fits 3e4 mask1.fits
    dmimgcalc infile=mask.fits infile2=mask1.fits op=mul clobber=yes outfile=mask.fits
    dmstat mask.fits centroid=no
done
