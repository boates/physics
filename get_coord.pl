#!/usr/bin/perl

# usage: get_coord.pl input rc
# where input is a file in the cnn format 
# and rc is a cut off radius within which number of neighbours 
# for each atom is counted 

$inputfile = shift(@ARGV);
$rc = shift(@ARGV);

$datafile = "coordination_rc_${rc}.dat";
$datafile1 = "av_coordination_rc_${rc}.dat";

open(IN, "<$inputfile");
chop(@Lines = <IN>);
close(IN);

$nh = 10;                   # number header lines
@header = @Lines[0 ... $nh-1];

foreach(@header) {
    chomp;
    if (/a =\s+(\S+)/) {$ca = $1;}      
    if (/b =\s+(\S+)/) {$cb = $1; }      # cell sizes
    if (/c =\s+(\S+)/) {$cc = $1; }
    if (/number_of_particles =\s+(\S+)/)   {$NA = $1; } # number_atoms
    if (/number_of_neighbours =\s+(\S+)/) {$NN = $1; }  # number neigh. listed
}

$NS = (@Lines - 10)/$NA;
for ($i=0;$i<= $NN;$i++) {$crd[$i] = 0;} # a list to store number of 0, 1, 2.. coordin. atoms

open(OUT,">$datafile");  open(OUT1,">$datafile1");
foreach(@header) {printf(OUT "$_ \n"); printf(OUT1 "$_ \n"); }

for ($is = 0; $is<$NS; $is++) {
    $i1 = $nh + $NA*$is; 
    $i2 = $i1 + $NA - 1;
    @snapshot = @Lines[$i1 ... $i2];   # isolate a snapshot

    
    foreach(@snapshot) {               # analyze atom by atom
	chomp;
	@el = (split /\s+/);           # get the data for the give atom        	
	$atom=$el[0]; $x0=$el[1]; $y0=$el[2]; $z0=$el[3];
	$atcoord = 0;

	$inn=0; $flag=1;

	while ($inn < $NN & $flag) {
	    $i = $inn + 4;
	    $natom = $el[$i];         # find which is the (inn+1)'th nearest neighbour
	    @el_n  = (split /\s+/, $snapshot[$natom]);
	    $x1 = $el_n[1]; $y1 = $el_n[2]; $z1 = $el_n[3]; # coordinates of that neigh.

	    $dx = $x1-$x0; $dx = $dx - $ca*&anint($dx/$ca);
	    $dy = $y1-$y0; $dy = $dy - $cb*&anint($dy/$cb);
	    $dz = $z1-$z0; $dz = $dz - $cc*&anint($dz/$cc);
	    
	    $r = sqrt($dx*$dx + $dy*$dy +$dz*$dz);  # distance to the neighbour
	    
	    if ($r < $rc) { $atcoord++} else { $flag = 0; }
	    $inn++;
	}
#	printf (OUT "$atcoord @el \n");   #use this to print all data in the cnn file
	printf (OUT "$atcoord \n");   #use this to print only the individual coordinations

	$crd[$atcoord]++; 
		       }
}

close(OUT);

for ($i=0;$i<= $NN;$i++) {$crd[$i] = $crd[$i]/$NS;}  # normalize

printf (OUT1 "average number of 0, 1, 2, ...coordinated atoms \n @crd \n");
close(OUT1);

sub anint {
    my $x = $_[0]; 
    my $i = int($x);
    if (abs($x-$i) >= 0.5) { 
        if ($x > 0){ $i = $i + 1; }
        if ($x < 0){ $i = $i - 1; }
    }
    $i;
}
