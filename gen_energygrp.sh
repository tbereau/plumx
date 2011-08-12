#!/bin/sh

# generate energygrps and energygrp_table variables for .mdp file.
# Append script output to .mdp file.
# do not include SC - SC (they use the generic table.xvg).
# read in 1 letter-code sequence to determine list of interactions.
# 
# name of amino acids. 
# Special 1-letter codes: 
#  - Z corresponds to end cap (CAP);
#  - B: Arginine (+); 
#  - J: Aspartic acid (-);
#  - O: Glutamic acid (-);
#  - U: Lysine (+);


die() {
    echo $@
    exit 1
}

[ -z $1 ] && die "Missing sequence file"

seq=`cat $1`
seq_length=`cat $1 | wc -c`

list_lipid=(CH PH GL ES AS AD AE N C)
list3aa=(ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM GLY HIS ILE LEU LYS0 LYSP MET PHE PRO SER THR TRP TYR VAL CAP)
list1aa=(A   R    B    N   D    J    C   Q   E    O    G   H   I   L   K    U    M   F   P   S   T   W   Y   V   Z  )
seq_aa=()

for (( i=0; i<$seq_length; ++i )); do
    letter=${seq:$i:1}
    if [ "$letter" != "" ]; then
	seq_index=-1
	for (( j=0; j<${#list1aa[@]}; ++j )); do
	    [ "$letter" == "${list1aa[$j]}" ] && seq_index=$j
	done
	[ "$seq_index" == "-1" ] && die "Amino acid code $letter is unrecognized."
	
	# if the residue hasn't been added to the list yet, add it now
	res_in_list=0
	for (( k=0; k<${#seq_aa[@]}; ++k )); do
	    [ "${seq_aa[$k]}" == "$seq_index" ] && res_in_list=1
	done
	[ "$res_in_list" == "0" ] && seq_aa=( ${seq_aa[@]-} $seq_index )
    fi
done

n_list_lipid=${#list_lipid[@]}
n_laa=${#list3aa[@]}

echo -ne "energygrps               = "
for ((i=0;i<$n_list_lipid;++i)); do
    echo -ne "${list_lipid[$i]} ";
done;
for (( k=0; k<${#seq_aa[@]}; ++k )); do
    echo -ne "${list3aa[${seq_aa[$k]}]} ";
done;


echo -ne "\nenergygrp_table          = "
for ((i=0;i<$n_list_lipid;++i)); do 
    for ((j=$i;j<$n_list_lipid;++j)); do 
	echo -ne "${list_lipid[$i]} ${list_lipid[$j]} "; 
    done; 
done; 

for ((i=0;i<$n_list_lipid;++i)); do 
    for (( k=0; k<${#seq_aa[@]}; ++k )); do
	echo -ne "${list3aa[${seq_aa[$k]}]} ${list_lipid[$i]} "; 
    done; 
done; 

echo ""
