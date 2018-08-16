#!/bin/sh
#
# Generate 'Timg_Gsyn_${freq}.fits' for each given frequency,
# using tool 'mkimg_galsyn_freq.py'.
#
# Aaron LI
# 2015/03/28
#

# Raw input radiation image at 408 MHz used as template.
RAWFILE="gal_sync_408_hi.fits"
# Spectral index map
INDEXFILE="gal_sync_index_hi.fits"

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
    OUTFILE="Timg_Gsyn_${f}.fits"
    mkimg_galsyn_freq.py --freq=${f} --rawfile="${RAWFILE}" \
        --indexfile="${INDEXFILE}" --outfile="${OUTFILE}" \
        --verbose --clobber --single
    echo ""
done

