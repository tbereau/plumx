#!/bin/sh
#
# Convert atomistic PDB file to CG PDB file.  Generates VMD script to convert
# PDB to PLUM-compatible input structure.
#
# CG model from
#       Bereau and Deserno, J. Chem. Phys. 130 (2009). 
#
# Tristan Bereau, 05/29/2012.
#

[ -z $1 ] && echo "Missing argument: file.pdb" && exit 1

pdb=$1
pdb_base=`basename $pdb .pdb`

cat > convert2cg.tmp.tcl <<EOF

mol load pdb $pdb

autopsf -mol 0

mol delete 0

if { ![file exists ${pdb_base}_autopsf.psf] } {
    puts "Error. autopsf produced an error."
    exit
}

mol load psf ${pdb_base}_autopsf.psf pdb ${pdb_base}_autopsf.pdb

set cg [atomselect top "not (sidechain or type HB HC or name OT2 OT3) or name CB"]

\$cg writepdb ${pdb_base}_cg.pdb
EOF

vmd < convert2cg.tmp.tcl >convert2cg.vmd.out

[ ! -f ${pdb_base}_cg.pdb ] && echo \
    "VMD failed to output CG structure. See convert2cg.vmd.out" \
    && exit 1

cat ${pdb_base}_cg.pdb 

rm ${pdb_base}_cg.pdb ${pdb_base}_autopsf.psf ${pdb_base}_autopsf.pdb \
    ${pdb_base}_autopsf.log convert2cg.tmp.tcl convert2cg.vmd.out
