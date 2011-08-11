#!/bin/sh
#
# Convert atomistic pdb file to CG pdb file. 
# Trim file and keep only atom names: N HN CA CB C O.
#
# CG model from 
#     Bereau and Deserno, J. Chem. Phys. 130 (2009).  
#
# Tristan Bereau, 09/08/2011.
#

die() {
    echo $@
    exit 1
}

[ -z $1 ] && die "Not enough arguments: 1) pdb file"

echo "REMARK  Generated with pdb_aa2cg.sh"
cat $1 | awk 'BEGIN{}
{
    if ($1=="ATOM" && ($3=="N" || $3=="HN" || $3=="CA" || $3=="CB" || $3=="C" || $3=="O")) {
        print $0;
    }
}'

