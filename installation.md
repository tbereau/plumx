# Installation #

## Gromacs ##


Please note that the following steps are NOT necessary if you want to use the
Lipid force field alone. The special kernel is needed for the Peptide and the
Protein-Lipid part of the force field.

## in case of Gromacs 4.6.6 ##

**Caution: Gromacs 4.6 patch does not run correctly. Please use Gromacs 4.5 (see below)**

  * Replace the charge-group--charge-group kernel to add an implementation of the hydrogen-bond interaction. To do so one can either apply the `nb_generic_cg-GMX-4.6.6.patch` to your gromacs source:
```
cd path/to/gromacs
patch -p1 < nb_generic_cg-GMX-4.6.6.patch
```
> or replace the charge-group--charge-group kernel in `src/gmxlib/nonbonded/nb_generic_cg.c` by the file included in this archive:
```
cp /path/to/plumx/nb_generic_cg-GMX-4.6.6.c path/to/gromacs/src/gmxlib/nonbonded/nb_generic_cg.c
```

  * recompile Gromacs
```
make
make install
```
## in case of Gromacs 4.5.5 ##

  * Replace the charge-group--charge-group kernel to add an implementation of the hydrogen-bond interaction. To do so one can either apply the `nb_generic_cg-GMX-4.5.4_5.patch` to your gromacs source:
```
cd path/to/gromacs
patch -p1 < nb_generic_cg-GMX-4.5.4_5.patch
```
> or replace the charge-group--charge-group kernel in `src/gmxlib/nonbonded/nb_generic_cg.c` by the file included in this archive:
```
cp /path/to/plumx/nb_generic_cg-GMX-4.5.4_5.c path/to/gromacs/src/gmxlib/nonbonded/nb_generic_cg.c
```

  * recompile Gromacs
```
make
make install
```

## in case of Gromacs 4.5.4 ##

  * Same as for Gromacs 4.5.5, except that: In the Gromacs directory, open `src/mdlib/ns.c`.  Around line 322, comment out the line:
```
gmx_fatal(FARGS,"The charge-group - charge-group force loops only \
  support systems with all intra-cg interactions excluded and no inter-cg \
  exclusions, this is not the case for this system.");
```
> This problem has been solved in Gromacs 4.5.5 and above.

  * recompile Gromacs
```
make
make install
```

Other versions of Gromacs do not yet have ready-made patches.

## PLUMX ##

Download PLUMX by following the instructions in https://code.google.com/p/plumx/source/checkout. You can either
  1. (easy solution) get a local copy of the code, allowing you updates (read-only)
```
git clone https://code.google.com/p/plumx/
```
  1. create your own clone, and thereby adapt the repository with an independent git repository. Note that any change will affect _your_ repository, not `plumx` directly. If you wish to contribute to `plumx`, request a code review.

The repository uses `git` to keep the code up to date. To download PLUMX, copy-paste the command shown to you at the aforementioned page in your terminal. This will create a subdirectory `plumx` and download a copy of the repository. To periodically take advantage of updates, you can simply run
```
git pull
```
to update the code to its latest version. Make sure, however, not to make any changes in the `plumx` subdirectories. Any simulation should be run in a separate directory (otherwise updates of the code will be problematic).

  * Export environment variable GMX\_NBLISTCG (see file 'source.me') to activate charge-group--charge-group kernel (this includes the Hbond interaction):
```
source source.me
```

> Note: in case of lipid-only simulation (i.e., no protein), it is best _not_ to export GMX\_NBLISTCG (or, alternatively, to set as 0) which is only required when the charge-group--charge-group kernel needs to use the C code (i.e., where the Hbond interaction is implemented).  This will speed up performance.