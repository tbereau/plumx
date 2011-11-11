#!/bin/sh
#
# Generate list of C6 C12 gromacs coefficients for PLUM protein-lipid force
# field (cross interactions only)
# Tristan Bereau (bereau@alumni.cmu.edu)
# Aug. 8, 2011


# list of side chains - last one is end cap (1-letter code: Z)
sidechains=(ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM GLY HIS ILE LEU LYS0
    LYSP MET PHE PRO SER THR TRP TYR VAL CAP)
# number of side chains
nr_sc=25
# side-chain max ID
sc_max=`expr "$nr_sc - 1" | bc -l`

# lists of sigma parameters for wca interactions. "-1" means lj rather than wca.
# (i.e., not wca). sigmas are in angstroms. convert to nm below.
# wca interactions with ch
sc_CH=(-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 4.0)
# wca interactions with ph
sc_PH=(9.54 9.32 6.66 7.45 7.42 9.27 8.32 7.74 7.69 9.94 8.73 9.28 9.92 10.25
    9.02 6.63 8.96 9.82 8.99 8.24 9.28 7.69 8.15 10.08 4.0)
# wca interactions with gl
sc_GL=(-1 -1 -1 -1 -1 -1 -1 -1 -1 5.71 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 4.0)
# wca interactions with e1, e2
sc_E1=(4.61 1.94 -1 3.01 1.79 2.39 -1 2.50 2.48 -1 4.34 0.60 6.73 6.73 -1 2.57
    -1 6.37 0.60 1.73 0.60 -1 1.98 1.87 4.0)
sc_E2=( ${sc_E1[@]} )
# wca interactions with as
sc_AS=(-1 4.17 4.17 3.86 2.74 -1 1.09 4.60 4.00 5.71 -1 4.15 -1 -1 1.48 4.74
1.47 -1 4.68 3.16 4.15 1.89 3.65 -1 4.0)
# wca interactions with ad
sc_AD=( ${sc_AS[@]} )
# wca interactions with ae
sc_AE=(-1 11.03 13.71 9.92 4.93 13.7 1.37 8.05 9.14 10.28 -1 7.74 -1 -1 1.19
10.67 1.18 -1 8.25 10.01 7.74 3.15 4.87 -1 4.0)

# lists of sigma parameters for lj interactions. "-1" means wca rather than
# lj. sigmas are in angstroms. convert to nm.
sclj_CH=(5.96 6.66 6.66 6.21 6.18 6.18 6.16 6.45 6.41 6.41 5.63 6.23 6.61 6.61
    6.63 6.63 6.59 6.77 6.20 5.97 6.23 6.99 6.79 6.42 -1)
sclj_PH=(-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
    -1 -1)
sclj_GL=(5.26 5.13 4.17 3.86 3.84 5.48 3.82 4.03 4.00 -1 4.44 4.98 5.91 5.91
    4.39 4.15 4.36 5.46 4.95 4.22 4.98 4.15 4.45 5.15 -1)
sclj_E1=(-1 -1 4.52 -1 -1 -1 0. -1 -1 0. -1 -1 -1 -1 0. -1 0. -1 -1 -1 -1
    0. -1 -1 -1)
sclj_E2=( ${sclj_E1[@]} )
sclj_AS=(4.21 -1 -1 -1 -1 0. -1 -1 -1 -1 3.94 -1 3.55 3.55 -1 -1 -1 3.64 -1 -1
    -1 -1 -1 3.43 -1) 
sclj_AD=( ${sclj_AS[@]} )
sclj_AE=(4.21 -1 -1 -1 -1 -1 -1 -1 -1 -1 3.94 -1 3.55 3.55 -1 -1 -1 3.64 -1 -1
    -1 -1 -1 3.43 -1)

# list of epsilon parameters for lj interactions. "-1" means wca rather than
# lj. epsilons are in kT, converto kJ/mol.
epslj_CH=(1.0 1.0 1.2 1.1 1.2 1.2 1.0 1.2 1.2 1.2 1.0 1.0 1.2 1.2 1.0 1.2 1.0
    1.2 1.0 1.0 1.0 1.2 1.2 1.2 -1)
epslj_PH=(-1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1 -1
-1 -1)
epslj_GL=(4.5 7.5 6.5 7.5 6.5 4.0 7.0 7.5 7.5 -1 4.5 7.5 5.5 5.5 7.5 6.5 7.5
    7.5 7.5 5.5 7.5 7.0 6.5 7.0 -1)
epslj_E1=(-1 -1 2.0 -1 -1 -1 0.0 -1 -1 0.0 -1 -1 -1 -1 0.0 -1 0.0 -1 -1 -1 -1
    0.0 -1 -1 -1)
epslj_E2=( ${epslj_E1[@]} )
epslj_AS=(0.85 -1 -1 -1 -1 0.0 -1 -1 -1 -1 0.7 -1 1.45 1.2 -1 -1 -1 1.1 -1 -1
    -1 -1 -1 1.0 -1)
epslj_AD=( ${epslj_AS[@]} )
epslj_AE=(0.85 -1 -1 -1 -1 -1 -1 -1 -1 -1 0.7 -1 1.5 1.2 -1 -1 -1 1.2 -1 -1 -1
    -1 -1 1.0 -1)


for i in CH PH GL E1 E2 AS AD AE; do
    scwca=`echo sc_${i}`
    eval scwca=( \${$scwca[@]} )
    sclj=`echo sclj_${i}`
    eval sclj=( \${$sclj[@]} )
    epslj=`echo epslj_${i}`
    eval epslj=( \${$epslj[@]} )

    for j in `seq 1 ${#sidechains[@]}`; do 
	# WCA by default - 2.58 factor corresponds to kT -> kJ/mol conversion
	c6=`expr "4*0.02*2.58" | bc -l`	    
	c12=0.
	if [ ${scwca[$j-1]} == "-1" ]; then
	    # LJ rather than WCA
	    c6=`expr "4. * 2.58 * ${epslj[$j-1]} * (${sclj[$j-1]}/10.)^6" | bc -l`
	    c12=`expr "4. * 2.58 * ${epslj[$j-1]} * (${sclj[$j-1]}/10.)^(12)" | bc -l`	   
	fi
	printf " $i   %4s   1     %6.5e   %6.5e\n" ${sidechains[$j-1]} $c6 $c12
    done
    
    # backbone
    for j in HBN N HBC; do 
	c6=`expr "4*0.02*2.58" | bc -l`	    
	c12=0.
	printf " $i   %4s   1     %6.5e   %6.5e\n" $j $c6 $c12
    done
done

