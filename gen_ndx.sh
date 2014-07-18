#!/bin/bash
#
# Generate index file for PLUM simulation.  
# Reads .gro file and runs make_ndx
# with a specific set of rules. This will create 
# groups compatible with the energy groups generated
# by gen_energygrp.sh.
#

die () {
    echo $1
    exit 1
}

[ -z $1 ] && die "Please provide .gro file"


file="tmp_ndx.XXX"
rm -f $file

nam_incr=1
# keep 0
echo "keep 0" >> $file
query=("a C CA O" "a N HN")
name=(C N)
# add side chains (except GLY, no CB)
for i in ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM HIS ILE LEU LYS0 \
    LYSP MET PHE PRO SER THR TRP TYR VAL CAP; do
    grep $i $1 > /dev/null
    if [ $? -eq 0 ]; then
	      query=( "${query[@]}" "a CB & r $i")
	      name=( "${name[@]}" $i)
    fi
done
query=("${query[@]}" "a CH" "a PH" "a GL" "a ES1 ES2" \
    "a AS11 AS12 AS13 AS14 AS21 AS22 AS24 AS25" "a AD23" "a AE15 AE26")
name=( "${name[@]}" CH PH GL ES AS AD AE)

n_query=${#query[@]}

for (( i=0; i<$n_query; ++i )); do
    echo "${query[$i]}" >> $file
    # [ "$i" -gt "0" ] && 
    echo "name $nam_incr ${name[$i]}" >> $file
    let nam_incr+=1
done


# for ((i=1;i<100;++i)); do
#     echo "del $i"
# done
echo "q" >> $file

cat $file
cat $file | make_ndx -f $1
rm -f $file
