#!/bin/sh
#
# Run some testings of cluster simulation.
# The generated cluster images will be used in the following
# power spectrum analysis for comparison.
#
# Aaron LI
# 2015/04/25
#

FREQ="050.00"
GPS_DAT="gps.dat"
TOOLS_DIR="../tools"

TESTS="$@"

DIR_INIT=`pwd -P`
for t in ${TESTS}; do
    T=`python -c "print('%02d' % int(${t}))"`
    DIR="test_${T}"
    echo "==========================================================="
    echo "*** ${DIR} ***"
    [ ! -d ${DIR} ] && mkdir ${DIR} && cd ${DIR}
    echo "Simulating gps data ..."
    #${TOOLS_DIR}/simulate_gps gps.dat
    ${TOOLS_DIR}/simulate_gps_trng gps.dat
    echo "Creating gps maps of frequency ${FREQ} MHz ..."
    ${TOOLS_DIR}/create_gps_map gps.dat ${FREQ} gps_${FREQ}.fits
    cd ${DIR_INIT}
done

