README - PointSrc/
==================

Weitian LI
Created: 2015/03/05
Updated: 2015/03/06


PROBLEMS
--------
1. `query_results/*`: where are them from? how to reproduce them?
2. CSS (compact steep spectrum) AGNs not considered here? Only GPS AGN?
3. `ptr${freq}.fits` & `gps_${freq}.fits`: differences? relations?


Data
----
* `query_results/*.query.result`:
  ???
  See above problems.


Tools
-----
* `gen_all_ptrs.sh`:
  1. Generate `gps.dat`, by invoking `simulate_gps.cc` to simulate
     the GPS (GHz peaked spectrum) AGN data. Then create maps of
     GPS AGN (`gps_${freq}.fits`) at given frequencies using the
     above generated `gps.dat`, by invoking `create_gps_map.cc`.
  2. Create maps of point source (`ptr${freq}.fits`), by
     invoking `make_map.cc`.
* `make_map.cc`:
  Generate maps of point source (`ptr*`) at given frequencies,
  providing the list of 'query result' files.
* `simulate_gps.cc`:
  Simulate the information of 'GPS (GHz peaked spectrum) AGN',
  which will be used in `create_gps_map.cc`.
* `create_gps_map.cc`:
  Create the maps of GPS AGNs using the above simulated information.

