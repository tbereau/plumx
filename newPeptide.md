# Work with different peptide(s) #

The following builds on the [tutorial example](helloWorld.md) and addresses specifically what to do to simulate a different peptide (whether with or without a membrane).

  * Write protein sequence file protein.seq: 1 protein per line, all 1-letter amino acid codes, no spaces between letters.  Example: 15-residue polyalanine:
```
   AAAAAAAAAAAAAAA
```
> Alternatively, one can extract the 1-letter sequence from a PDB file using:
```
   ./extract.1aa.seq.sh protein.pdb
```
> note however that amino acids with different protonation states will need to be written by hand (replace "?" signs).

> Special amino acid codes:
    1. Z: end cap (N-terminal or C-terminal group)
    1. B: Arginine [+]; (neutral arginine is R)
    1. J: Aspartic acid [-]; (neutral aspartic acid is D)
    1. O: Glutamic acid [-]; (neutral glutamic acid is E)
    1. U: Lysine [+]; (neutral lysine is K)

  * Generate protein force field file by using script:
```
./plum_prot_gen_itp.pl protein.seq > prot.itp
```

  * Write topology file topol.top which combines the lipid, protein, and plum `.itp` files.  (See example file: `popc72_walp/topol.top`).

  * Generate `.gro` file from a pdb structure which contains the lipid and protein coordinates (concatenated) from the script
```
   ./pdb2gro.sh initial.pdb Lx Ly Lz
```
> > where initial.pdb is the pdb structure, and Lx, Ly, and Lz are the three box sizes (in nm) in the x, y, and z directions, respectively.


> Note 1: an atomistic protein can simply be coarse-grained by using the script
```
./pdb_aa2cg.sh protein.pdb
```

> Note 2: pdb2gmx does _not_ work with PLUM because the residue database has not yet been implemented (see .rtp files in the Gromacs manual). Pay attention to nonconventional names for some residues (e.g., histidine should be HIS).

  * Copy grompp.mdp file from example directory popc72\_walp. Use the script
```
./gen_energygrp.sh [-no_lipid] protein.seq 
```
to generate 'energygrps' and 'energygrp\_table' variables. Use the [-no\_lipid] option in case you're running a protein simulation without the membrane.  Insert script output in grompp.mdp file.  Make sure the 'userint{1,2,3}' variables are set as in the popc72\_walp example (they are required for proper Hbond calculation).  Also, do not change the order of [atomtypes](atomtypes.md) in `plum.itp` file (the userint{2,3} variables correspond to the atomtypes of beads HBN and HBC, hydrogen-bonding capable amide and carbonyl groups).  The optional -no\_lipid option removes lipid parameters in case none are present in the simulation box.

  * Check consistency of .gro file: all amino acid codes should be one of the following:
```
   ALA ARG0 ARGP ASN ASP0 ASPM CYS GLN GLU0 GLUM GLY HIS ILE LEU LYS0
    LYSP MET PHE PRO SER THR TRP TYR VAL CAP
```
> > where a trailing 0 corresponds to the neutral state of the amino acid, while a trailing "P" and "M" correspond to PLUS and MINUS charges, respectively. CAP corresponds to the capping group for transmembrane proteins. A small sed script can easily replace, say, ASP for ASP0 amino acids (if all ASP residues are to be neutralized in the simulation):
```
   sed 's/ASP /ASP0/g' conf.gro
```

  * Generate index file index.ndx from .gro file:
```
./gen_ndx.sh conf.gro
```