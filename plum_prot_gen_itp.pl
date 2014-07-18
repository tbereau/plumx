#!/usr/bin/perl
#
# Generate Gromacs topology file (.itp) for PLUM peptide model from given
# amino acid sequence (1-letter code string).
#
# Tristan Bereau (08/02/2011)
#

# name of amino acids. 
# Special 1-letter codes: 
#  - Z corresponds to end cap;
#  - B: Arginine (+); 
#  - J: Aspartic acid (-);
#  - O: Glutamic acid (-);
#  - U: Lysine (+);
our %aa_name = (
    G=>"GLY",  A=>"ALA",  D=>"ASP0", N=>"ASN", E=>"GLU0", Q=>"GLN", V=>"VAL",
    L=>"LEU",  I=>"ILE",  M=>"MET",  T=>"THR", S=>"SER",  C=>"CYS", K=>"LYS0",
    R=>"ARG0", H=>"HIS",  F=>"PHE",  P=>"PRO", W=>"TRP",  Y=>"TYR", Z=>"CAP",
    B=>"ARGP", J=>"ASPM", O=>"GLUM", U=>"LYSP");

our %aa_mass = (
    G=>" 75.0", A=>" 89.0", D=>"133.0", N=>"132.0", E=>"147.0", Q=>"146.0", V=>"117.0",
    L=>"131.0", I=>"131.0", M=>"149.0", T=>"119.0", S=>"105.0", C=>"121.0", K=>"146.0",
    R=>"174.0", H=>"155.0", F=>"165.0", P=>"115.0", W=>"204.0", Y=>"181.0", Z=>" 75.0",
    B=>"174.0", J=>"133.0", O=>"147.0", U=>"146.0");

$prog_name="plum_prot_gen_itp.pl";

# parse arguments
$numargs = scalar @ARGV;
if ($numargs != 1) {
    print STDERR "Sequence file must be specified.\n";
    exit;
} else {
    $file_seq = <$ARGV[0]>;
}
$VERBOSE=1;

# read sequence(s)
open SEQ, '<', $file_seq or die "Can't open sequence file: $!";
@sequences = ();
$seq1 = "";
@begin_chain = ();
@end_chain = ();
$number_of_chains = 0;
while ($line = <SEQ>) {
    chomp $line;
    $seq1 = $line;
    push @sequences, $seq1;
    $number_of_chains++;
}

# print header
print "; topology file for CG PLUM protein simulation.\n";
print "; Generated from $prog_name\n";
print ";\n; ** Amino acid sequence **\n";
for (my $i=0; $i < $number_of_chains; $i++) 
{
    print "; Chain ",$i+1,":";
    print " $sequences[$i]\n";
}

# print moleculetype
print "[ moleculetype ]\n";
print "; molname    nrexcl\n";
print "  Protein    2\n\n";

# print [ atoms ]
my $bead_count=1;
my $res_count=1;
my @type = ("HBN", "  H", " CA", "   ", "HBC", "  O");
my $type_pro = "N";
my @name = (" N", "HN", "CA", "CB", " C", " O");
my @name_pro = (" N", "CA", "CB", " C", " O");
my @name_gly = (" N", "HN", "CA", " C", " O");
my @cgnr = (1,1,2,3,4,4);
my $cgnr_count=0;
my @mass = (14.0,1.0,12.0,0.0,12.0,16.0);
print "[ atoms ]\n";
print ";    nr   type    resnr   residue atom    cgnr charge      mass\n";
for (my $i=0; $i < $number_of_chains; $i++)
{
    my @seqi = split(//, $sequences[$i]);
    for (my $j=0; $j < @seqi; $j++)
    {
        for (my $k=0; $k < 6; $k++)
        {
            if ($k == 0 && $seqi[$j] eq "P") {
                printf(" %6d    %3s     %4d       %3s   %2s   %5d %6.3f    %6.2f\n",
                       $bead_count,$type_pro,
                       $res_count,PRO,$name[$k],
                       $cgnr[$k]+$cgnr_count,0,$mass[$k]);		
            } elsif ($k == 1 && $seqi[$j] eq "P") {
                # Do nothing
                $bead_count--;
            } elsif ($k != 3) {
                printf(" %6d    %3s     %4d       %3s   %2s   %5d %6.3f    %6.2f\n",
                       $bead_count,$type[$k],
                       $res_count,$aa_name{$seqi[$j]},$name[$k],
                       $cgnr[$k]+$cgnr_count,0,$mass[$k]);
            } else {
                if ($seqi[$j] eq "G") {
                    # Do nothing
                    $bead_count--;
                } else {
                    printf(" %6d    %3s     %4d       %3s   %2s   %5d %6.3f    %6.2f\n",
                           $bead_count,$aa_name{$seqi[$j]},
                           $res_count,$aa_name{$seqi[$j]},$name[$k],
                           $cgnr[$k]+$cgnr_count,0,$aa_mass{$seqi[$j]});
                }
            }
            $bead_count++;
        }
        $res_count++;
        $cgnr_count+=4;
    }
}

# print [ bonds ]
$res_count=0;
$at_count=0;
my @ai = (1,1,3,3,5,5);
my @aj = (2,3,4,5,6,7);
my @aj2 = (2,3,4,5,6,1);
my @ai_gly  = (1,1,0,3,4,4);
my @ai2_gly = (1,1,0,3,5,5);
my @aj_gly  = (2,3,0,4,5,6);
my @aj2_gly = (2,3,0,5,6,1);
my @ai_pro  = (0,1,2,2,4,4);
my @ai2_pro = (0,1,3,3,5,5);
my @aj_pro  = (0,2,3,4,5,6);
my @aj2_pro = (0,3,4,5,6,1);
my $kbond=77400; # in units of kJ/mol - see CG model paper.
# bond list: NH NCa CaCb CaC' C'O (C'N)
my @bondlist = (0.1, 0.1455, 0.1530, 0.1510, 0.1235, 0.1325);
print "\n[ bonds ]\n";
print ";    ai        aj    funct     b0             kb\n";
for (my $i=0; $i < $number_of_chains; $i++)
{
    my @seqi = split(//, $sequences[$i]);
    for (my $j=0; $j < @seqi; $j++)
    {
        my $max_k=5;
        # is it the end of the chain?
        if ($j < @seqi-1) { 
            $max_k = 6;	
        }
        for (my $k=0; $k < $max_k; $k++)
        {
            if (($seqi[$j] eq "G" and $k == 2) or ($seqi[$j] eq "P" and $k == 0)) {
                # Do nothing
            } elsif ($seqi[$j] eq "G") {
                printf(" %6d    %6d     %4d ;   %6.5e    %6.5e   %2s - %2s interaction\n",
                       $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],1,@bondlist[$k],$kbond,@name[@ai2_gly[$k]-1],@name[@aj2_gly[$k]-1]);
            } elsif ($seqi[$j] eq "P") {
                printf(" %6d    %6d     %4d ;   %6.5e    %6.5e   %2s - %2s interaction\n",
                       $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],1,@bondlist[$k],$kbond,@name[@ai2_pro[$k]-1],@name[@aj2_pro[$k]-1]);
            } else {
                printf(" %6d    %6d     %4d ;   %6.5e    %6.5e   %2s - %2s interaction\n",
                       $at_count+$ai[$k],$at_count+$aj[$k],1,@bondlist[$k],$kbond,@name[@ai[$k]-1],@name[@aj2[$k]-1]);
            }
        }
        $res_count++;
        if ($seqi[$j] eq "G" or $seqi[$j] eq "P") {
            $at_count += 5;
        } else {
            $at_count += 6;      
        }
    }
}

$res_count=0;
$at_count=0;
# print [ angles ]
print "\n[ angles ]\n";
print ";    ai        aj         ak    funct     k0             ka\n";
my @ai = (2,1,1,4,3,3,5);
my @aj = (1,3,3,3,5,5,7);
my @ak = (3,4,5,5,6,7,9);
my @aj2 = (1,3,3,3,5,5,1);
my @ak2 = (3,4,5,5,6,1,3);
my @ai_gly  = (2,0,1,0,3,3,4);
my @aj_gly  = (1,0,3,0,4,4,6);
my @ak_gly  = (3,0,4,0,5,6,8);
my @ak_glypro = (3,0,4,0,5,6,7);
my @ai_pro  = (0,1,1,3,2,2,4);
my @aj_pro  = (0,2,2,2,4,4,6);
my @ak_pro  = (0,3,4,4,5,6,8);
my @ak_xxxpro = (3,4,5,5,6,7,8);
my @ak_propro = (0,3,4,4,5,6,7);
my $kangle = 774.0; # in units of kJ/mol
my @anglelist = (115,108,111,113,122,116,122);
for (my $i=0; $i < $number_of_chains; $i++)
{
    my @seqi = split(//, $sequences[$i]);
    for (my $j=0; $j < @seqi; $j++)
    {
        my $max_k=5;
        # is it the end of the chain?
        if ($j < @seqi-1) {
            $max_k = 7;
        }
        for (my $k=0; $k < $max_k; $k++)
        {
            if (($seqi[$j] eq "G" and ($k == 1 or $k==3)) or ($seqi[$j] eq "P" and $k == 0)) {
                # Do nothing
            } elsif ($k == 6 and $seqi[$j] eq "G" and $seqi[$j+1] eq "P") {
                printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                       $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_glypro[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
            } elsif ($k == 6 and $seqi[$j] eq "P" and $seqi[$j+1] eq "P") {
                printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                       $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_propro[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
            } elsif ($k == 6 and $seqi[$j+1] eq "P") {
                printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                       $at_count+$ai[$k],$at_count+$aj[$k],$at_count+$ak_xxxpro[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
            } else {
                if ($seqi[$j] eq "G") {
                    printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                           $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_gly[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
                } elsif ($seqi[$j] eq "P") {
                    printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                           $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_pro[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
                } else {
                    printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                           $at_count+$ai[$k],$at_count+$aj[$k],$at_count+$ak[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
                }
            }
        }
        
        $res_count++;
        if ($seqi[$j] eq "G" or $seqi[$j] eq "P") {
            $at_count += 5;
        } else {
            $at_count += 6;      
        }
    }
}

$res_count=0;
$at_count=0;
# print [ dihedrals ]
print "\n[ dihedrals ]\n";
print ";    ai        aj         ak         al    funct   phi        cp            mult\n";
# N 1; H 2; CA 3; CB 4; C 5; O 6
my @ai         = (1,3,1, 5,5,3);
my @aj         = (3,5,3, 7,7,5);
my @ak         = (5,7,5, 9,8,6);
my @al         = (4,9,7,11,9,7);
my @aj2        = (3,5,3, 1,1,5);
my @ak2        = (5,1,5, 3,2,6);
my @al2        = (4,3,1, 5,3,1);
my @ai_gly     = (0,3,1, 4,4,3);
my @aj_gly     = (0,4,3, 6,6,4);
my @ak_glyxxx  = (0,6,4, 8,7,5);
my @ak_glypro  = (0,6,4, 7,0,5);
my @al_xxxgly  = (4,9,7,10,9,7);
my @al_glygly  = (0,8,6, 9,8,6);
my @al_glyxxx  = (0,8,6,10,8,6);
my @al_glypro  = (0,7,6,10,0,6);
my @ai_pro     = (1,2,1, 4,4,2);
my @aj_pro     = (2,4,2, 6,6,4);
my @ak_proxxx  = (4,6,4, 8,7,5);
my @ak_xxxpro  = (5,7,5, 8,0,6);
my @ak_propro  = (4,6,4, 7,0,5);
my @al_proxxx  = (3,8,6,10,8,6);
my @al_xxxpro  = (4,8,7,10,8,7);
my @al_progly  = (3,8,6, 9,8,6);
my @al_propro  = (3,7,6, 9,7,6);
# omega dihedral=172.86. All cases except proline.
my @kdih = (43.86,172.86,-0.774,0.774,150.,150.);
my $komega_pro = 7.74;
my @dihedrallist = (60,0,180,180,0.,0.);
my @multiplicities = (1,1,1,1,1,1);
my $mult_omega_pro = 2;
for (my $i=0; $i < $number_of_chains; $i++)
{
    my @seqi = split(//, $sequences[$i]);
    for (my $j=0; $j < @seqi; $j++)
    {
        my $max_k=1;
        # is it the end of the chain?
        if ($j < @seqi-1) {
            $max_k = 6;
        }
        for (my $k=0; $k < $max_k; $k++)
        {
            if ($seqi[$j] eq "G") {
                if ($k != 0) {
                    if ($seqi[$j+1] eq "P") {
                        if ($k != 4) {
                            printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                                   $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_glypro[$k],$at_count+$al_glypro[$k],1,
                                   @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                        }
                    } elsif ($seqi[$j+1] eq "G") {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_glyxxx[$k],$at_count+$al_glygly[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    } else {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_gly[$k],$at_count+$aj_gly[$k],$at_count+$ak_glyxxx[$k],$at_count+$al_glyxxx[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    }
                }
            } elsif ($seqi[$j] eq "P") {
                if ( $k==1) {
                    # omega of proline is specific
                    if ($seqi[$j+1] eq "P") {
                        if ($k != 4 ) {
                            printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                                   $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_propro[$k],$at_count+$al_propro[$k],1,
                                   @dihedrallist[$k],$komega_pro,$mult_omega_pro,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                        }
                    } else {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_proxxx[$k],$at_count+$al_proxxx[$k],1,
                               @dihedrallist[$k],$komega_pro,$mult_omega_pro,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);                        
                    }
                } else {
                    if ($seqi[$j+1] eq "P") {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_propro[$k],$at_count+$al_propro[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    } elsif ($seqi[$j+1] eq "G") {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_proxxx[$k],$at_count+$al_progly[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    } else {
                        printf(" %6d    %6d     %6d     %6d     %4d ;   %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai_pro[$k],$at_count+$aj_pro[$k],$at_count+$ak_proxxx[$k],$at_count+$al_proxxx[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    }
                }
            } else {
                if ($seqi[$j+1] eq "G") {
                    printf(" %6d    %6d     %6d     %6d     %4d ;  %6.2f    %+6.5e    %2d;   %2s - %2s - %2s - %2s interaction\n",
                           $at_count+$ai[$k],$at_count+$aj[$k],$at_count+$ak[$k],$at_count+$al_xxxgly[$k],1,
                           @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);   
                } elsif ($seqi[$j+1] eq "P") {
                    if ($k != 4) {
                        printf(" %6d    %6d     %6d     %6d     %4d ;  %6.2f    %+6.5e    %2d;   %2s - %2s - %2s - %2s interaction\n",
                               $at_count+$ai[$k],$at_count+$aj[$k],$at_count+$ak_xxxpro[$k],$at_count+$al_xxxpro[$k],1,
                               @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
                    }
                } else {
                    printf(" %6d    %6d     %6d     %6d     %4d ;  %6.2f    %+6.5e    %2d;   %2s - %2s - %2s - %2s interaction\n",
                           $at_count+$ai[$k],$at_count+$aj[$k],$at_count+$ak[$k],$at_count+$al[$k],1,
                           @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);                    
                }
            }
        }
        $res_count++;
        if ($seqi[$j] eq "G" or $seqi[$j] eq "P") {
            $at_count += 5;
        } else {
            $at_count += 6;      
        }
    }
}

