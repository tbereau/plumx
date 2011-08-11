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
	    } elsif ($k != 3) {
		printf(" %6d    %3s     %4d       %3s   %2s   %5d %6.3f    %6.2f\n",
		       $bead_count,$type[$k],
		       $res_count,$aa_name{$seqi[$j]},$name[$k],
		       $cgnr[$k]+$cgnr_count,0,$mass[$k]);
	    } else {
		printf(" %6d    %3s     %4d       %3s   %2s   %5d %6.3f    %6.2f\n",
		       $bead_count,$aa_name{$seqi[$j]},
		       $res_count,$aa_name{$seqi[$j]},$name[$k],
		       $cgnr[$k]+$cgnr_count,0,$aa_mass{$seqi[$j]});
	    }
	    $bead_count++;
	}
	$res_count++;
	$cgnr_count+=4;
    }
}

# print [ bonds ]
$res_count=0;
my @ai = (1,1,3,3,5,5);
my @aj = (2,3,4,5,6,7);
my @aj2 = (2,3,4,5,6,1);
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
	    printf(" %6d    %6d     %4d ;   %6.5e    %6.5e   %2s - %2s interaction\n",
		   $res_count*6+$ai[$k],$res_count*6+$aj[$k],1,@bondlist[$k],$kbond,@name[@ai[$k]-1],@name[@aj2[$k]-1]);
	}
	$res_count++;
    }
}

$res_count=0;
# print [ angles ]
print "\n[ angles ]\n";
print ";    ai        aj         ak    funct     k0             ka\n";
my @ai = (2,1,1,4,3,3,5);
my @aj = (1,3,3,3,5,5,7);
my @ak = (3,4,5,5,6,7,9);
my @aj2 = (1,3,3,3,5,5,1);
my @ak2 = (3,4,5,5,6,1,3);
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
            printf(" %6d    %6d     %6d     %4d ;   %6.5e    %6.5e   %2s - %2s - %2s interaction\n",
                   $res_count*6+$ai[$k],$res_count*6+$aj[$k],$res_count*6+$ak[$k],1,@anglelist[$k],$kangle,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1]);
        }
        $res_count++;
    }
}

$res_count=0;
# print [ dihedrals ]
print "\n[ dihedrals ]\n";
print ";    ai        aj         ak         al    funct   phi        cp            mult\n";
# N 1; H 2; CA 3; CB 4; C 5; O 6
my @ai = (1,3,1,5,5,3);
my @aj = (3,5,3,7,7,5);
my @ak = (5,7,5,9,8,6);
my @al = (4,9,7,11,9,7);
my @aj2 = (3,5,3,1,1,5);
my @ak2 = (5,1,5,3,2,6);
my @al2 = (4,3,1,5,3,1);
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
	    # omega of proline is specific
	    if ( $k==1 && $seqi[$j] eq "P") {
		printf(" %6d    %6d     %6d     %6d     %4d    %6.2f    %+6.5e    %2d;  %2s - %2s - %2s - %2s interaction\n",
		       $res_count*6+$ai[$k],$res_count*6+$aj[$k],$res_count*6+$ak[$k],$res_count*6+$al[$k],1,
		       @dihedrallist[$k],$komega_pro,$mult_omega_pro,@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
	    } else {
		printf(" %6d    %6d     %6d     %6d     %4d ;  %6.2f    %+6.5e    %2d   %2s - %2s - %2s - %2s interaction\n",
		       $res_count*6+$ai[$k],$res_count*6+$aj[$k],$res_count*6+$ak[$k],$res_count*6+$al[$k],1,
		       @dihedrallist[$k],$kdih[$k],$multiplicities[$k],@name[@ai[$k]-1],@name[@aj2[$k]-1],@name[@ak2[$k]-1],@name[@al2[$k]-1]);
	    }
        }
        $res_count++;
    }
}

$res_count=0;
# print [ exclusions ]
print "\n[ exclusions ]\n";
print "; explicitly exclude all nonbonded HN-O interactions that are 3 bonds apart.\n";
print ";    ai        aj\n";
for (my $i=0; $i < $number_of_chains; $i++)
{
    my @seqi = split(//, $sequences[$i]);
    for (my $j=0; $j < @seqi; $j++)
    {
        # is it the end of the chain?
        if ($j < @seqi-1) {
	    printf(" %6d    %6d ;  %2s - %2s interaction\n",
		   $res_count*6+6,$res_count*6+8,"O","H");
        }
        $res_count++;
    }
}


