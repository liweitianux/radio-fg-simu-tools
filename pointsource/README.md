Simulation of point sources
===========================


Dependencies
------------
* `fio`: see `../fio`
* `CCfits`: CFITSIO C++ bindings
  (`apt install libccfits-dev`)


Build
-----
* `make`


Tools
-----
* `make_ptr_map`:
  Generate maps of point sources at given frequencies based on the
  simulation results released by [Wilman et al. 2008](http://adsabs.harvard.edu/abs/2008MNRAS.388.1335W).
* `simulate_gps`:
  Simulate the properties of *GPS (GHz peaked spectrum)* and
  *CSS (Compact steep spectrum)* AGNs, which will be used by
  `create_gps_map` below.
* `create_gps_map`:
  Create the maps of GPS and CSS AGNs.
* `gen_all_ptrs.sh`:
  Script to help simulate all the point sources.


Data
----
The `make_ptr_map` tool requires the simulation data from
[Wilman et al. 2008](http://adsabs.harvard.edu/abs/2008MNRAS.388.1335W),
which are published at [The SKA Simulated Skies - SEX](http://s-cubed.physics.ox.ac.uk/s3_sex).

We obtained the simulation data of 100 deg^2 in 2010, and
the dataset is quite huge (2.4 GB after bzip2 compression).
We can provide a temporary download link upon request.


Usage
-----
An example script to simulate the ptr maps:

```sh
#!/bin/sh

case "$1" in
    ""|-h|--help)
        echo "usage: ${0##*/} <freqMHz> ..."
        exit 1
        ;;
esac

FOV=10  # [deg]
PIX_SIZE=20  # [arcsec]
IMG_SIZE=$(echo "${FOV} * 3600 / ${PIX_SIZE}" | bc)
PREFIX="ptr${IMG_SIZE}fov${FOV}_"
DBLIST="wilman2008_db.list"
OUTCSV="${PREFIX}sources.csv"

echo "FoV: ${FOV} [deg]"
echo "Pixel size: ${PIX_SIZE} [arcsec]"
echo "Image size: ${IMG_SIZE}"
echo "Prefix: ${PREFIX}"
echo "Frequencies: $@"

make_ptr_map \
    -f ${FOV} -s ${IMG_SIZE} \
    -i ${DBLIST} -O ${OUTCSV} -o ${PREFIX} "$@"
```
