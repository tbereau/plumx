#!/bin/sh
#
# Generate entire set of tables used for PLUM protein-lipid force field
# (cross interactions only)
#  Tristan Bereau (bereau@alumni.cmu.edu)
#  Aug. 5, 2011 

if [ -z $1 ]; then
    echo "Error missing argument: distance column. Usage:"
    echo "$0 dist.dat"
    exit
fi

# SC-SC interaction goes into the generic table: table.xvg
echo "Generating table SC - SC"
cat $1 | awk 'BEGIN{
    # side chain sigma
    s=0.50;
    # 2**(1/6)*sigma
    rc=0.56123
    # max value
    max=1e7;
    
    print "# PLUM force field: protein side chain - protein side chain interaction"
    print "# C6  column: multiply by 4*eps_hp, where eps_hp = 4.5 (see paper)"
    print "# C12 column: multiply by 4*eps_hp*eps_ij, where eps_ij is the normalized"
    print "# Lorentz-Berthelot mixing coefficient."
    print "#"
    print "# see paper: Bereau and Deserno, JCP 130, 235106 (2009)"
    print "#"
    print "# dist           U_Coulomb       F_Coulomb       U_vdwC6         F_vdwC6         U_vdwC12        F_vdwC12"
    
}{
    if ($1==0) 
	printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,max,0,-1./4,0;
    else if ($1<rc) { 
	lj=(s/$1)**12 - (s/$1)**6 + 0.25;
	if (lj > max) lj = max;    
	printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,lj,0,-1./4,0;
    } else 
	printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,0.,0,(s/$1)**12 - (s/$1)**6,0;
}' > table.xvg

die () {
    echo $1;
    exit 1;
}

# LJ table
gen_lj () {
    [ -z $1 ] && die "missing parameter"
    cat $1 | awk 'BEGIN{
        # max value
	max=1e7;
	
	print "# PLUM force field: 0, -1/r6, 1/r12 table"
	}{
	if ($1==0) 
	    printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,-1.*max,0,max,0;
	else { 
		c6 =(1/$1)**6;
		c12=(1/$1)**12;
		if (c6 > max ) c6 = max;
		    if (c12 > max ) c12 = max;
			printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,-1.*c6,0,c12,0;
            }
    }'    
}

# WCA table
gen_wca () {
    [ -z $2 ] && die "missing parameter"
    cat $1 | awk -v s=$2 'BEGIN{
    # sigma
    #s=$2;
    # 2**(1/6)*sigma
    rc=2.**(1/6)*s;
    # max value
    max=1e7;
	
    printf "# PLUM force field: WCA interaction with sigma = %s\n",s
    }
    {
	if ($1==0) 
	    printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,max,0,0,0;
	else if ($1<rc) { 
		    lj=(s/$1)**12 - (s/$1)**6 + 0.25;
		    if (lj > max) lj = max;    
		        printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,lj,0,0,0;
		} else 
	        printf "%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\t%6.5e\n", $1,0,0,0,0,0,0;
    }'
}

sN=0.29;
sC=0.36;
sSC=0.50;

# Interactions with N
for i in N C SC; do 
    sX=`echo s${i}`
    eval sX=\$$sX
    s=`echo "0.5 * (${sN} + ${sX})" | bc -l`
    echo "Generating table  N - $i"
    gen_wca $1 $s > table_N_${i}.xvg
done

# Interactions with C
for i in C SC; do 
    sX=`echo s${i}`
    eval sX=\$$sX
    s=`echo "0.5 * (${sC} + ${sX})" | bc -l`
    echo "Generating table  C - $i"
    gen_wca $1 $s > table_C_${i}.xvg
done

# list of side chains - last one is end cap (1-letter code: Z)
sidechains=(ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM GLY HIS ILE LEU LYS0
    LYSP MET PHE PRO SER THR TRP TYR VAL CAP)
# number of side chains
nr_sc=25
# side-chain max ID
sc_max=`expr "$nr_sc - 1" | bc -l`

# lists of sigma parameters for WCA interactions. "-1" means LJ rather than WCA.
# (i.e., not WCA). Sigmas are in Angstroms. Convert to nm below.
# WCA interactions with CH
sc_CH=(-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 4.0)
# WCA interactions with PH
sc_PH=(9.54 9.32 6.66 7.45 7.42 9.27 8.32 7.74 7.69 9.94 8.73 9.28 9.92 10.25
    9.02 6.63 8.96 9.82 8.99 8.24 9.28 7.69 8.15 10.08 4.0)
# WCA interactions with GL
sc_GL=(-1 -1 -1 -1 -1 -1 -1 -1 -1 5.71 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 4.0)
# WCA interactions with E1, E2
sc_ES=(4.61 1.94 -1 3.01 1.79 2.39 -1 2.50 2.48 -1 4.34 0.60 6.73 6.73 -1 2.57
    -1 6.37 0.60 1.73 0.60 -1 1.98 1.87 4.0)
# WCA interactions with AS
sc_AS=(-1 4.17 4.17 3.86 2.74 -1 1.09 4.60 4.00 5.71 -1 4.15 -1 -1 1.48 4.74
1.47 -1 4.68 3.16 4.15 1.89 3.65 -1 4.0)
# WCA interactions with AD
sc_AD=( ${sc_AS[@]} )
# WCA interactions with AE
sc_AE=(-1 11.03 13.71 9.92 4.93 13.7 1.37 8.05 9.14 10.28 -1 7.74 -1 -1 1.19
10.67 1.18 -1 8.25 10.01 7.74 3.15 4.87 -1 4.0)

for i in CH PH GL ES AS AD AE; do
    list=`echo sc_${i}`
    eval list=( \${$list[@]} )
    # loop over sidechains
    for j in `seq 0 $sc_max`; do 
	echo "Generating table $i - ${sidechains[$j]}"
	if [ "${list[$j]}" == "-1" ]; then
	    gen_lj  $1 > table_${i}_${sidechains[$j]}.xvg
	else
	    gen_wca $1 `echo "${list[$j]}/10." | bc -l` > table_${i}_${sidechains[$j]}.xvg
	fi
    done
    # loop over backbone beads
    for j in N C; do
	sX=`echo s${j}`
	eval sX=\$$sX
	s=`echo "0.5 * ${sX} + 0.3" | bc -l`
	echo "Generating table $i - $j"
	gen_wca $1 $s > table_${i}_${j}.xvg
    done
done

for i in ${sidechains[@]}; do
    for j in N C; do
	sX=`echo s${j}`
        eval sX=\$$sX
        s=`echo "0.5 * ($sSC + ${sX})" | bc -l`
        echo "Generating table $i - $j"
        gen_wca $1 $s > table_${j}_${i}.xvg
    done
done

echo "Generated all tables into files table*.xvg"
