#!/usr/bin/perl -w
use strict;
use warnings;

#Three files are needed: 1) sin from Wang et al. 2) 3did_flat file 3) list of variations
if ($#ARGV != 2) {die "Program used with parameters [SIN] [3did] [list of variations]\n";}

#	query interface from ipfam or 3did and from pisa or ???
#	print structuraly resolved interactions
#	interpret interfaces
#		pfam, functional residues, complexation residues		

#from list of variations
	#map variations on the network
	#identify pairs or chains of  interactions = > group 
#from list of proteins
#	get interactions
#	get domain pairs (of interacting protein pair)
#	get list of interacting AA

#parse_variations
my @variations=open_file($ARGV[2]);
my (@variat, %seen, @uniq);
for(my $i=0;$i<$#variations+1;$i++){
	my @tab = split(/\t/,$variations[$i]);#Gene_symbol, Entrez, Ref AA, Position, Var AA, Ensemble_transcript, Uniprot AC, Iniprot ID
	
	push (@{$variat[$i]}, @tab);
    	$seen{$tab[1]}++;
#print "$variat[$i][1]\n";
}	
@uniq = keys %seen;

#get interactions from Wang (structural and/or PPI)
my @wang_sin=open_file($ARGV[0]);
my (@sin);
#my @wang_pin=open_file($ARGV[2]);	
for(my $i=0;$i<$#wang_sin+1;$i++){#parse interactions
	my @tab=split(/\t/,$wang_sin[$i]);#protA protB pfamA seq_str_A seq_end_A pfamB seq_str_B seq_end_B
	push (@{$sin[$i]}, @tab);
#print "kaka $sin[$i][0]\n";
}


#get interfaces from 3did
my @threedid=open_file($ARGV[1]);
my $str="NO";
my $ni=0;
my (@int,@interactom,$nc,@cmap);
for(my $i=0;$i<$#threedid+1;$i++){#parse interactions
	my @int = split(/\t/,$threedid[$i]);
	my $id = $int[0];
	$id =~ s/\W//;	
	if($id eq "ID"){
		$ni++;
		$interactom[$ni][0]=$int[1];#name domain 1
		$interactom[$ni][1]=$int[2];#name domain 2
		$interactom[$ni][2]=substr($int[3],0,index($int[3],"."));#PFAM 1 id	
		$interactom[$ni][3]=substr($int[4],0,index($int[4],"."));#PFAM 1 id
		$str="NO";
	} 
	elsif($id eq "3D"){
		$nc=0;#reset contact map counter
		$interactom[$ni][4]=$int[1];#pdb id
		my @tab1=split(":",$int[2]);
		$interactom[$ni][5]=$tab1[0];#domain 1 chain id
		my @tab2=split("-",$tab1[1]);
		$interactom[$ni][6]=$tab2[0];#domain 1 start 
		$interactom[$ni][7]=$tab2[1];#domain 1 end
		my @tab3=split(":",$int[3]);
		$interactom[$ni][8]=$tab3[0];#domain 2 chain id
		my @tab4=split("-",$tab3[1]);
		$interactom[$ni][9]=$tab4[0];#domain 2 start
		$interactom[$ni][10]=$tab4[1];#domain 2 end
		$interactom[$ni][11]=$int[4];#score
		$interactom[$ni][12]=$int[5];#zscore
		$str="YES";		
	}	
	else{
	 if($str eq "YES"){
		$nc++;#count interacting AA pairs in given interface
		$cmap[$ni][$nc][0]=$int[0];#AA1 id
		$cmap[$ni][$nc][1]=$int[1];#AA2 id	
		$cmap[$ni][$nc][2]=$int[2];#AA1
		$cmap[$ni][$nc][3]=$int[3];#AA2	
		$cmap[$ni][$nc][4]=$int[4];#interaction type (side chain or main chain)
	 }
	}
}

#map variations
my(%subsy,@waryjaty,@din,%dvar,%cout);
my $k=0;
for (my $j=0;$j<$#wang_sin+1;$j++){
	my $dom1=$sin[$j][0]."_".$sin[$j][2];
	my $dom2=$sin[$j][1]."_".$sin[$j][5];
	$dvar{$dom1}=0;
	$dvar{$dom2}=0;
	my $inter=$dom1."\t".$dom2;	
	$cout{$dom1}=0;
	$cout{$dom2}=0;
	for(my $i=0;$i<$#variat+1;$i++){
		if($variat[$i][1] eq $sin[$j][0]){
		  $k++;
		  my @tab=@{$sin[$j]};
		  my @tab2=@{$variat[$i]};
    		  $subsy{$inter}++;
		  #$waryjaty[$i][6]=$sin[$j][2];
		  #przypisywac id "reakcji"?
		  push (@{$din[$k]}, @tab);
		  my $pos = $variat[$i][3];
		  my $ds =$sin[$j][3];
		  my $de =$sin[$j][4];	
		  if($pos ge $ds and $pos le $de){
			$dvar{$dom1}++;
		  }		
		  else{$cout{$dom1}++;}
#print "kaka $din[$k][0] $din[$k][1] $din[$k][2] $din[$k][3] $din[$k][4] $din[$k][5] $din[$k][6] $din[$k][7]\n";
		}
		elsif($variat[$i][1] eq $sin[$j][1]){
		#elsif($uniq[$i] eq $sin[$j][1]){
		  $k++;
		  my @tab=@{$sin[$j]};
    		  $subsy{$inter}++;
		  push (@{$din[$k]}, @tab);
		  my $pos = $variat[$i][3];
		  my $ds =$sin[$j][6];
		  my $de =$sin[$j][7];
		  if($pos ge $ds and $pos le $de){
			$dvar{$dom2}++;
		  }		
		  else{$cout{$dom2}++;}
		}
	#print #"waraiacja /pfam /pfam domain interakting/ interfejs/ functionaly important group (part of pfam domain) /klasyfikacja (w domenie, w interfejsie,  
	}
}

my $filec = substr($ARGV[2],0,rindex($ARGV[2],"."));
my $file = $filec."-graph.txt";
open(OUT, "> $file") or die "Can not open an output file: $!";
foreach my $klucz (keys %subsy){
	my $r1=0;
	my $r2=0;
	printf OUT "$klucz\n";	
	my @tab=split(/\s+/,$klucz);
	my $t1=$dvar{$tab[0]}+$cout{$tab[0]};#total number of variations 
	my $t2=$dvar{$tab[1]}+$cout{$tab[1]};
	$r1=($dvar{$tab[0]}/$t1) unless $t1 == 0;
	$r2=($dvar{$tab[1]}/$t2) unless $t2 == 0;
	print "$tab[0]\t$dvar{$tab[0]}\t$cout{$tab[0]}\t$t1\t$r1\t$tab[1]\t$dvar{$tab[1]}\t$cout{$tab[1]}\t$t2\t$r2\n";
}

print "$k\n";


#map variations on interfaces
#@subsy=split(/\s+/,$variat[$i][4]);

sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
		
