#!/bin/sh
#
# Generate 'Timg_Gff_${freq}.fits' for each given frequency,
# using tool 'mkimg_galff_freq.py'.
#
# Aaron LI
# 2015/03/28
#

# Input Halpha map used as the template for Galactic free-free radiation map.
INFILE="halpha_1024.fits"

errmsg() {
    echo "$@" > /dev/stderr
}

if [ $# -ne 1 ]; then
    errmsg "Usage: `basename $0` <freq_list>"
    exit 1
fi

freq_list="$1"
cat "${freq_list}" | while read f; do
    echo "freq: ${f}"
    OUTFILE="Timg_Gff_${f}.fits"
    mkimg_galff_freq.py --freq=${f} --infile="${INFILE}" \
        --outfile="${OUTFILE}" --verbose --clobber
    echo ""
done

