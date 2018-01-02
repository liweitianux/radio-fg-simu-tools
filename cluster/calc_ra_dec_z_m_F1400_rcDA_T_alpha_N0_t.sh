#!/bin/sh
#
# Calculate the total cluster number by randomly selecting
# data from "dndmdata" with 'calc_totalnum'.
#
# Aaron LI
# 2015/04/12
#

if [ $# -ne 2 ] && [ $# -ne 3 ]; then
    printf "Usage:\n"
    printf "    `basename $0` MassLowerLimit B_Local_uG [clobber_yN]\n"
    exit 1
fi

# mass lower limit of clusters (total mass)
# default: 2.0e14 (Msun)
MASS_LOWER_LIMIT="$1"
# B_local: local magnetic field intensity (unit: uG)
B_LOCAL="$2"

# clobber
case "$3" in
    [Yy]*|[Tt]*)
        CLOBBER="YES"
        ;;
    *)
        CLOBBER="NO"
        ;;
esac
printf "CLOBBER=${CLOBBER}\n"

TOOLS_DIR="./tools"
ZMF_FILE="dndmdata"

# dark matter ratio
DM_RATIO="0.83"
# lower limit of dark matter mass
DM_MASS_LOWER_LIMIT=`python -c "print(${DM_RATIO}*${MASS_LOWER_LIMIT})"`

printf "Calculating total number of clusters ...\n"
CMD="${TOOLS_DIR}/calc_totalnum ${DM_MASS_LOWER_LIMIT} ${ZMF_FILE}"
printf "CMD: ${CMD}\n"
OUTPUT=`${CMD}`
TOTAL_NUMBER=`echo "${OUTPUT}" | awk '{ print $2 }'`
printf "MASS_LOWER_LIMIT=${MASS_LOWER_LIMIT}\n"
printf "DM_MASS_LOWER_LIMIT=${DM_MASS_LOWER_LIMIT}\n"
printf "TOTAL_NUMBER=${TOTAL_NUMBER}\n"

printf "Simulating (RA, DEC) data ...\n"
RA_DEC_DAT="ra_dec.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${RA_DEC_DAT}" ]; then
    CMD="${TOOLS_DIR}/sphxy ${TOTAL_NUMBER} ${RA_DEC_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Simulating alpha data ...\n"
ALPHA_DAT="alpha.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${ALPHA_DAT}" ]; then
    CMD="${TOOLS_DIR}/alpha_sim ${TOTAL_NUMBER} ${ALPHA_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Simulating evolution time (age) data ...\n"
AGE_DAT="age.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${AGE_DAT}" ]; then
    CMD="${TOOLS_DIR}/age_sim.py ${TOTAL_NUMBER} ${AGE_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Simulating z & m data ...\n"
Z_M_DAT="z_m.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${Z_M_DAT}" ]; then
    CMD="${TOOLS_DIR}/z_m_simulation ${DM_MASS_LOWER_LIMIT} ${TOTAL_NUMBER} ${ZMF_FILE} ${Z_M_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Calculating T & P1400 data ...\n"
Z_M_T_P1400_DAT="z_m_T_P1400.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${Z_M_T_P1400_DAT}" ]; then
    CMD="${TOOLS_DIR}/calc_T_P1400 ${Z_M_DAT} ${Z_M_T_P1400_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Calculating Lrosat data ...\n"
Z_M_T_P1400_LROSAT_DAT="z_m_T_P1400_Lrosat.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${Z_M_T_P1400_LROSAT_DAT}" ]; then
    CMD="${TOOLS_DIR}/f_rosat/calc_Lrosat ${Z_M_T_P1400_DAT} ${Z_M_T_P1400_LROSAT_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Calculating F1400 & rcDA data ...\n"
Z_M_F1400_RCDA_T_DAT="z_m_F1400_rcDA_T.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${Z_M_F1400_RCDA_T_DAT}" ]; then
    CMD="${TOOLS_DIR}/calc_F1400_rcDA ${Z_M_T_P1400_LROSAT_DAT} ${Z_M_F1400_RCDA_T_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Pasting data ...\n"
Z_M_F1400_RCDA_T_ALPHA_DAT="z_m_F1400_rcDA_T_alpha.dat"
RA_DEC_Z_M_F1400_RCDA_T_ALPHA_DAT="ra_dec_z_m_F1400_rcDA_T_alpha.dat"
paste -d" " ${Z_M_F1400_RCDA_T_DAT} ${ALPHA_DAT} > ${Z_M_F1400_RCDA_T_ALPHA_DAT}
paste -d" " ${RA_DEC_DAT} ${Z_M_F1400_RCDA_T_ALPHA_DAT} > ${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_DAT}

printf "Calculating synchrotron emission coefficient (NN0) data ...\n"
RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_DAT="ra_dec_z_m_F1400_rcDA_T_alpha_NN0.dat"
if [ "${CLOBBER}" = "YES" ] || [ ! -f "${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_DAT}" ]; then
    CMD="${TOOLS_DIR}/calc_NN0 ${B_LOCAL} ${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_DAT} ${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_DAT}"
    printf "CMD: ${CMD}\n"
    ${CMD}
else
    printf "Skipped!\n"
fi

printf "Pasting data (${AGE_DAT}) ...\n"
RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_AGE_DAT="ra_dec_z_m_F1400_rcDA_T_alpha_NN0_age.dat"
paste -d" " ${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_DAT} ${AGE_DAT} > ${RA_DEC_Z_M_F1400_RCDA_T_ALPHA_NN0_AGE_DAT}

