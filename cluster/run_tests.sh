#!/bin/sh
#
# Run some testings of cluster simulation.
# The generated cluster images will be used in the following
# power spectrum analysis for comparison.
#
# Aaron LI
# 2015/04/25
#

B_LOCAL=2
MASS_LOWER="2e14"
FREQ="050.00"

TESTS="$@"

DIR_INIT=`pwd -P`
for t in ${TESTS}; do
    T=`python -c "print('%02d' % int(${t}))"`
    DIR="test_${T}"
    echo "==========================================================="
    echo "*** ${DIR} ***"
    [ ! -d ${DIR} ] && mkdir ${DIR} && cd ${DIR}
    ln -s ../dndmdata ../calc_*.sh ../tools .
    echo "${FREQ}" > f.list
    ./calc_ra_dec_z_m_F1400_rcDA_T_alpha_N0_t.sh ${MASS_LOWER} ${B_LOCAL}
    ./calc_Fnu_mkimg.sh ${B_LOCAL} ra_dec_z_m_F1400_rcDA_T_alpha_NN0_age.dat f.list
    cd ${DIR_INIT}
done

