# Hello Peptide Simulation #

Here's a quick guide to run your first peptide simulation.

  * Copy the example directory `/PATH/TO/plumx/decaalanine` out of your local repository.  The directory contains everything required to run a `(Ala)_{10}` peptide.


  * Copy `plum_tables/*.xvg` into simulation directory.

  * The 1-letter protein sequence has been saved in `protein.seq`.

  * Copy `plum.itp` to simulation directory.

  * Generate protein force field file by using script:
```
./plum_prot_gen_itp.pl protein.seq > prot.itp
```

  * The topology file has been saved in `topol.top`.

  * The initial configuration file has been saved in `conf.gro`.

  * The gromacs input script has been saved in `grompp.mdp`.

  * The index file has been saved in `index.ndx`.

  * Generate the gromacs input files
```
grompp -f grompp.mdp -c init.gro -p topol.top -n index.ndx -norenum
```
Note that the `-norenum` option is **necessary** to ensure that the hydrogen-bond potential works correctly.

  * **Make sure your GMX\_NBLISTCG environment variable is set!**

  * Simulate
```
mdrun -v
```