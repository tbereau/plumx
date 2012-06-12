#!/bin/sh
#
# Convert atomistic PDB file to CG PDB file.  Generates VMD script to convert
# PDB to PLUM-compatible input structure.
#
# CG model from
#       Bereau and Deserno, J. Chem. Phys. 130 (2009). 
#
# Tristan Bereau, 05/29/2012.

[ -z $1 ] && echo "Missing argument: file.pdb" && exit 1

pdb=$1
pdb_base=`basename $pdb .pdb`

cat > convert2cg.tmp.tcl <<EOF

mol load pdb $pdb

autopsf -mol 0

mol delete 0

if { ![file exists ${pdb_base}_autopsf.psf] } {
    puts "Error. autopsf produced an error."
    exit
}

mol load psf ${pdb_base}_autopsf.psf pdb ${pdb_base}_autopsf.pdb

set cg [atomselect top "not (sidechain or type HB HC or name OT2 OT3) or name CB HT1"]
set firsth [atomselect top "name HT1"]
\$firsth set name HN
set lasto [atomselect top "name OT1"]
\$lasto set name O

\$cg writepdb ${pdb_base}_cg.pdb
EOF

vmd < convert2cg.tmp.tcl >convert2cg.vmd.out

[ ! -f ${pdb_base}_cg.pdb ] && echo \
    "VMD failed to output CG structure. See convert2cg.vmd.out" \
    && exit 1

# Reorder atoms
cat > reorder.tmp.pl <<EOF
#!/usr/bin/env perl
#
# Reorder atoms in PDB
#

@at_nm = ("N ","HN","CA","CB","C ","O ");
\$chain = "P1";


\$id_j = 0;

@store_lines = ();
\$last_i = 1;

while (defined(\$line = <STDIN>)) {
  if (\$line  =~ /HETATM/ || \$line =~ /ATOM/) {
    \$at_id_i = substr \$line, 4, 7;
    \$at_id_n = substr \$line, 13,2;
    \$at_id_res = substr \$line, 17,3;
    \$at_chain = substr \$line, 72,2;
    if (\$at_chain !~ \$chain) {
      # If we change of chain, look for stored residues
      for (\$ele = 0; \$ele < @store_lines; \$ele++) {
        \$ele_i = substr \$store_lines[\$ele], 4, 7;
        \$ele_n = substr \$store_lines[\$ele], 13, 2;
        \$ele_chain = substr \$store_lines[\$ele], 72, 2;
        if (\$ele_n =~ \$at_nm[\$id_j % 6] && \$ele_chain =~ \$at_chain) {
          \$at_id_res = substr \$store_lines[\$ele],17,3; 
          \$f_last_i = sprintf("%7d",\$last_i);
          \$newline = substr \$store_lines[\$ele], 4, 7, \$f_last_i;
          print "\$store_lines[\$ele]";
          \$at_id_res = substr \$store_lines[\$ele],17,3;
          \$chain = substr \$store_lines[\$ele], 72, 2;
          delete \$store_lines[\$ele];
          \$last_i += 1;
          \$id_j += 1;
          # Special cases for GLY and PRO
          if (\$at_id_res =~ "GLY" && \$id_j % 6 == 3) {
            \$id_j += 1;
          }
          if (\$at_id_res =~ "PRO" && \$id_j % 6 == 1) {
            \$id_j += 1;
          }      
        }
      }
    }

    if (\$at_nm[\$id_j % 6] !~ \$at_id_n) {
      for (\$ele = 0; \$ele < @store_lines; \$ele++) {
        \$ele_n = substr \$store_lines[\$ele], 13, 2;
        \$ele_chain = substr \$store_lines[\$ele], 72, 2;
        if (\$ele_n =~ \$at_nm[\$id_j % 6] && \$ele_chain =~ \$chain) {
          \$f_last_i = sprintf("%7d",\$last_i);
          \$newline = substr \$store_lines[\$ele], 4, 7, \$f_last_i;
          print "\$store_lines[\$ele]";
          \$at_id_res = substr \$store_lines[\$ele],17,3;
          \$chain = substr \$store_lines[\$ele], 72, 2;
          delete \$store_lines[\$ele];
          \$last_i += 1;
          \$id_j += 1;
          # Special cases for GLY and PRO
          if (\$at_id_res =~ "GLY" && \$id_j % 6 == 3) {
            \$id_j += 1;
          }
          if (\$at_id_res =~ "PRO" && \$id_j % 6 == 1) {
            \$id_j += 1;
          }
        }
      }
      push(@store_lines,\$line);
    } else {
      if (\$at_id_i !~ \$last_i) {
        \$f_last_i = sprintf("%7d",\$last_i);
        \$newline = substr \$line, 4, 7, \$f_last_i;
      }
      print \$line;
      \$last_i += 1;
      \$id_j += 1;
      # Special cases for GLY and PRO
      if (\$at_id_res =~ "GLY" && \$id_j % 6 == 3) {
        \$id_j += 1;
      }
      if (\$at_id_res =~ "PRO" && \$id_j % 6 == 1) {
        \$id_j += 1;
      }
    }
  } else {
    if (\$line !~ /END/) {
      print "\$line";
    }
  }
}

while (@store_lines > 0) {
  \$found_one = 0;
  for (\$ele = 0; \$ele < @store_lines; \$ele++) {
    \$ele_n = substr \$store_lines[\$ele], 13, 2;
    \$ele_chain = substr \$store_lines[\$ele], 72, 2;
    if (\$at_nm[\$id_j % 6] =~ \$ele_n) {
      \$f_last_i =sprintf("%7d",\$last_i);
      \$newline = substr \$store_lines[\$ele], 4, 7,\$f_last_i;
      print "\$store_lines[\$ele]";
      delete \$store_lines[\$ele];
      \$last_i += 1;
      \$id_j += 1;
      \$found_one = 1;
      # Special cases for GLY and PRO
      if (\$at_id_res =~ "GLY" && \$id_j % 6 == 3) {
        \$id_j += 1;
      }
      if (\$at_id_res =~ "PRO" && \$id_j % 6 == 1) {
        \$id_j += 1;
      }
    }
  }
  if (\$found_one == 0) {
    print "Error. Can't properly reorder the atoms. See PDB or contact\n";
    print "Tristan Bereau (bereau@alumni.cmu.edu)\n";
  }
}

print "END\n";

EOF

cat ${pdb_base}_cg.pdb | perl reorder.tmp.pl

#rm ${pdb_base}_cg.pdb ${pdb_base}_autopsf.psf ${pdb_base}_autopsf.pdb \
#    ${pdb_base}_autopsf.log convert2cg.tmp.tcl convert2cg.vmd.out reorder.tmp.pl
