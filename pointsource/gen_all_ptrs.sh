#!/bin/sh

TOOLS_DIR="./tools"

if [ $# -lt 1 ]; then
    echo "Usage:"
    echo "    $0 <freq_list_MHz>"
    exit 1
fi

if [ -f "gps.dat" ]; then
    echo "Using existing 'gps.dat'."
else
    echo "Simulating gps data ..."
    ${TOOLS_DIR}/simulate_gps gps.dat
fi

for i in $@; do
    echo "Creating gps maps of frequency $i MHz ..."
    ${TOOLS_DIR}/create_gps_map gps.dat $i gps_${i}.fits
done

echo "Creating ptr maps ..."
${TOOLS_DIR}/make_map query_results.list $@ ptr_

