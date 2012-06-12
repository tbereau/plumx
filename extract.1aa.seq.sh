#!/bin/bash
#
# Reads in pdb file; extract 1-letter-code amino-acid sequence.  
#
# Tristan Bereau (05/29/2012)

list3aa=(ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM GLY HIS ILE LEU LYS0 \
    LYSP MET PHE PRO SER THR TRP TYR VAL CAP)
list1aa=(A   R    B    N   D    J    C   Q   E    O    G   H   I   L   K    \
    U    M   F   P   S   T   W   Y   V   Z  )

die() {
    echo $@
    exit 1
}

[ -z $1 ] && die "Missing PDB file"

pdb=`cat $1 | awk '/CA/{if (\$0~"ATOM" || \$0~"HETATM") print substr(\$0,18,3)}'`
pdb_l=`echo \$pdb | wc -w`
let pdb_l*=4
echo "Reading 3-letter code sequence from $1:" >&2
echo $pdb >&2
echo "" >&2
echo "For charged residues, the following codes apply:" >&2
echo "  Z: end cap (CAP)" >&2
echo "  B: ARG+" >&2
echo "  J: ASP-" >&2
echo "  O: GLU-" >&2
echo "  U: LYS+" >&2
echo "Please make substitutions accordingly." >&2
echo "" >&2
echo "Final sequence (output to .seq file):" >&2

one_missing_res=0
for (( i=0; i<$pdb_l; i+= 4 )); do
    letter=${pdb:$i:3}
    found_res=0
    for (( j=0; j<${#list3aa[@]}; ++j )); do
        [ "$letter" == "${list3aa[$j]}" ] && echo -ne "${list1aa[$j]}" \
            && found_res=1
    done
    [ "$found_res" == "0" ] && echo -ne "?" && one_missing_res=1
    
done

echo -ne "\n"
[ "$one_missing_res" == "1" ] && \
    echo "Error: one or more amino acids were not recognized." >&2