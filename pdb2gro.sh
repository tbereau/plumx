#!/bin/sh
#
# Convert pdb file to gro file. Arguments:
#  - pdb file
#  - box size X (in nm)
#  - box size Y (in nm)
#  - box size Z (in nm)
#
#

die() {
    echo $@
    exit 1
}

[ -z $4 ] && die "Not enough arguments: 1) pdb file 2) box size X 3) box size Y 4) box size Z"

echo "Generated with pdb2gro.sh"
cat $1 | awk 'BEGIN{c=0}
{
    if ($1=="ATOM") {
	c++;
    }
}
END{printf "%5d\n",c;}'

cat $1 | awk 'BEGIN{res=1; restemp=1; npart=1}
{
    if ($1=="ATOM") {
        if ($6 != restemp) { res++; restemp = $6}
	printf "%5d%-4s  %4s%5d %7.3f %7.3f %7.3f\n", res, $4, $3, npart, $7/10., $8/10., $9/10.;
        npart++;
        if (npart == 100000) {
            npart = 0;
        }
        if (res == 100000) {
            res = 0;
        }
    }
}'
printf " %9.5f %9.5f %9.5f\n" $2 $3 $4;
