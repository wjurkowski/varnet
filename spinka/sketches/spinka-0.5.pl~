#!/usr/bin/perl -w
use strict;
use warnings;

#Three files are needed: 1-2) sin and pin from Wang et al. 3) list of variations
if ($#ARGV != 1) {die "Program used with parameters [list of variations] [SIN]\n";}

#parse_variations
my @variations=open_file($ARGV[0]);
my(@variat, %seen, @uniq);
for(my $i=0;$i<$#variations+1;$i++){
	my @tab = split(/\t/,$variations[$i]);#Gene_symbol, Entrez, UniprotAC, uniprotID, Ref AA, Position, Var AA
	push (@{$variat[$i]}, @tab);
    	$seen{$tab[1]}++;
}	
@uniq = keys %seen;

#get interactions from Wang sin (or other but in same format)
#entrez ID should already be mapped to gene symbols
my @wang_sin=open_file($ARGV[1]);
my (@sin);	
for(my $i=0;$i<$#wang_sin+1;$i++){#parse interactions
	my @tab=split(/\t/,$wang_sin[$i]);#protA pfamA seq_str_A seq_end_A protB pfamB seq_str_B seq_end_B
	push (@{$sin[$i]}, @tab);
}

#map variations and prepare DIN network
my(%din,%node1,%node2,%dvar,%cout);
my $k=0;
for (my $j=0;$j<$#wang_sin+1;$j++){
	my $dom1=$sin[$j][0]."_".$sin[$j][1];
	my $dom2=$sin[$j][4]."_".$sin[$j][5];
	my $inter=$dom1."\t".$dom2;	
	$dvar{$dom1}=0;#
	$dvar{$dom2}=0;
	$cout{$dom1}=0;
	$cout{$dom2}=0;
	for(my $i=0;$i<$#variat+1;$i++){
		if(($variat[$i][0] eq $sin[$j][0]) or ($variat[$i][0] eq $sin[$j][4])){
		  $node1{$inter}=$dom1;
		  $node2{$inter}=$dom2;	
		  $din{$inter}++;
		}
		if($variat[$i][0] eq $sin[$j][0]){
		  $k++;
		  my $pos = $variat[$i][5];
		  my $ds =$sin[$j][2];
		  my $de =$sin[$j][3];	
		  if($pos ge $ds and $pos le $de){
			$dvar{$dom1}++;#count variations in domain 1
		  }		
		  else{$cout{$dom1}++;}#count variations outside domain 1
#print "dupa $sin[$j][0] $node1{$inter} $dvar{$dom1} $cout{$dom1} \n";
		}
		elsif($variat[$i][0] eq $sin[$j][4]){
		  $k++;
		  my $pos = $variat[$i][5];
		  my $ds =$sin[$j][6];
		  my $de =$sin[$j][7];
		  if($pos ge $ds and $pos le $de){
			$dvar{$dom2}++;#count variations in domain 2
		  }		
		  else{$cout{$dom2}++;}#count variations outside domain 2
		} 
	}
}


my $filec = substr($ARGV[0],0,rindex($ARGV[0],"."));
my $file = $filec."-din-graph.txt";
my $file3 = $filec."-din-nodeAttr.txt";
my $file4 = $filec."-din-edgeAttr.txt";
open(OUT, "> $file") or die "Can not open an output file: $!";
open(OUT3, "> $file3") or die "Can not open an output file: $!";
open(OUT4, "> $file4") or die "Can not open an output file: $!";
printf OUT4 "edge\tedgeGV\ttotal_edge\tr_edge\n";#print edge attributes

my(%hasGV,%GVinD,%GVoutD,%nattri);
my $edgeGV=0;
foreach my $klucz (keys %din){
	my $r1=0;
	my $r2=0;
	my $rpair=0;

	printf OUT "$node1{$klucz}\t$node2{$klucz}\n";#print DIN network

	my @tab=split(/\s+/,$klucz);
	my $t1=$dvar{$tab[0]}+$cout{$tab[0]};#total number of variations 
	my $t2=$dvar{$tab[1]}+$cout{$tab[1]};#total number of variations
	$r1=($dvar{$tab[0]}/$t1) unless $t1 == 0;
	$r2=($dvar{$tab[1]}/$t2) unless $t2 == 0;
	my $tpair=$t1+$t2;
	my $dpairvar=$dvar{$tab[0]}+$dvar{$tab[1]};
	$rpair=($dpairvar/$tpair) unless $tpair == 0;

	#SNV inside interacting domains 
	if($dvar{$tab[0]} > 0){
		$GVinD{$tab[0]}=1;
	}
	else{
		$GVinD{$tab[0]}=0;
	}	
	if($dvar{$tab[1]} > 0){
		$GVinD{$tab[1]}=1;
	}
	else{
                $GVinD{$tab[1]}=0;
        }

	#SNV outside interacting domains
	if($cout{$tab[0]} > 0){
		$GVoutD{$tab[0]}=1;
	}
	else{
		$GVoutD{$tab[0]}=0;
	}
	if($cout{$tab[1]} > 0){
		$GVoutD{$tab[1]}=1;
	}
	else{
		$GVoutD{$tab[1]}=0;
	}
	

	if(($dvar{$tab[0]} > 0) or ($cout{$tab[0]} > 0)){
		$hasGV{$tab[0]}=1;
		$edgeGV=1;
	}
	else{$hasGV{$tab[0]}=0;}
	if(($dvar{$tab[1]} > 0) or ($cout{$tab[1]} > 0)){
		$hasGV{$tab[1]}=1;
		$edgeGV=1;
	}
	else{$hasGV{$tab[1]}=0;}
	if(($hasGV{$tab[0]} == 1) and ($hasGV{$tab[1]} == 1)){$edgeGV = 2;}	
	#printf OUT3 "$tab[0]\t$GVinD{$tab[0]}\t$GVoutD{$tab[0]}\t$t1\t$r1\n";#print attributes of first node
	#printf OUT3 "$tab[1]\t$GVinD{$tab[1]}\t$GVoutD{$tab[1]}\t$t1\t$r2\n";#print attributes of the second node
	$nattri{$tab[0]}="$node1{$klucz}\t$GVinD{$tab[0]}\t$GVoutD{$tab[0]}\t$t1\t$r1\n";#print attributes of first node
	$nattri{$tab[1]}="$node2{$klucz}\t$GVinD{$tab[1]}\t$GVoutD{$tab[1]}\t$t2\t$r2\n";#print attributes of the second node
	printf OUT4 "$node1{$klucz} (pp) $node2{$klucz}\t$edgeGV\t$tpair\t$rpair\n";#print edge attributes	
	#print some statistics
#	Domain 1;#SNV in domain;#SNV outside dom.;col2/total numb.;Domain 2;#SNV in domain;#SNV outside dom.;col6/total numb.;total on domains/total in edge    	
#	print "$tab[0]\t$dvar{$tab[0]}\t$cout{$tab[0]}\t$r1\t$tab[1]\t$dvar{$tab[1]}\t$cout{$tab[1]}\t$r2\t$rpair\n";
}
#print node attributes
printf OUT3 "node1\tGVinD\tGVoutD\tTotal\tr_GVinD\n";
foreach my $klucz (keys %nattri){printf OUT3 $nattri{$klucz};}
#print "$k\n";



sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file $file_name : $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
		