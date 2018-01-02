#!/bin/sh
#
# Calculate F${nu}, S0 data, and make image of cluster radiation component.
#
# Aaron LI
# 2015/04/22
#

TOOLS_DIR="./tools"

errmsg() {
    echo "$@" > /dev/stderr
}

if [ $# -ne 3 ]; then
    errmsg "Usage: `basename $0` <B_local_uG> <data: ra_dec_z_m_F1400_rcDA_T_alpha_NN0_age.dat> <freq_list>"
    exit 1
fi

B_LOCAL="$1"
DATAFILE="$2"
FREQ_LIST="$3"

cat "${FREQ_LIST}" | while read f; do
    echo "======================================="
    echo "freq: ${f}"
    echo "   calculating F${f} ..."
    OUTFILE="ra_dec_z_m_F${f}_rcDA_T_alpha.dat"
    ${TOOLS_DIR}/calc_Fnu_t ${B_LOCAL} ${f} ${DATAFILE} ${OUTFILE}
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "   calculating S0 & beta ..."
    INFILE="${OUTFILE}"
    OUTFILE="ra_dec_F${f}_S0_rcut_rcDA_beta.dat"
    ${TOOLS_DIR}/calc_s0 ${INFILE} ${OUTFILE}
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "   generating image ..."
    INFILE="${OUTFILE}"
    OUTFILE="img_cl_${f}.fits"
    ${TOOLS_DIR}/gen_2d/gen_2d_img_cut ${INFILE} ${OUTFILE}
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    echo "   converting to brightness image ..."
    INFILE="${OUTFILE}"
    OUTFILE="T${INFILE}"
    ${TOOLS_DIR}/transfer_sur ${f} ${INFILE} ${OUTFILE}
    echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
done

