#
# This XSPEC/TCL script is used to generate the conversion table between
# the bolometric flux (0.1-20 keV) and soft flux (0.5-2.0 keV, for ROSAT).
#
# The generated table 'fconv_table.inc' will be included into 'f_rosat.cc'.
#
# Usage:
#     $ xspec - gen_fconv_table.tcl
#
# Aaron LI
# 2015/04/12
#

dummyrsp 0.1 20.0
abund grsa
model mekal & 1 & 1e-6 & 0.3 & 0 & /*

set ff [ open fconv_table.inc w ]
for {set z 0} {$z < 7.5} {set z [expr $z+0.1]} {
    for {set T 0.5} { $T < 20 } {set T [expr $T+0.1]} {
        newpar 4 $z
        newpar 1 $T
        flux 0.1 20.0
        tclout flux
        scan $xspec_tclout "%f" bolo_flux
        flux 0.5 2.0
        tclout flux
        scan $xspec_tclout "%f" soft_flux
        puts $ff "$z,  $T,  $bolo_flux,    $soft_flux,"
    }
}

tclexit
