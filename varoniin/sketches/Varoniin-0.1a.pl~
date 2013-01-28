#!/usr/bin/perl -w
use strict;
use warnings;

#Four files are needed: 1) list of variations 2) sin and pin files from Wang et al. 3) interactome file
if ($#ARGV != 1) {die "Program used with parameters [list of variations] [interactome]\n";}

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
my @variations=open_file($ARGV[0]);
my(@variat, %seen, @uniq);
for(my $i=0;$i<$#variations+1;$i++){
	my @tab = split(/\t/,$variations[$i]);#Gene_symbol, Entrez, Ref AA, Position, Var AA, Ensemble_transcript, Uniprot AC, Iniprot ID
	
	push (@{$variat[$i]}, @tab);
    	$seen{$tab[1]}++;
#print "$variat[$i][1]\n";
}	
@uniq = keys %seen;

#for given AA position of a variant 
#	get corresponding pdb ID
#	check hash of given pdb

#read pfam - pdb mapping from pfam website: pfam definitions according to pdb resid numbering
#	for the moment we can assume that residue numbers in PDB are same as in uniprot seq.

#read in interactome
my @interactome=open_file($ARGV[1]);
my(@waryjaty,%cmapAAisGV,%varAAc,%intfc,%intfc_1,%intfc_2,%varAAc_ss,%varAAc_sm,%varAAc_mm,%varAAc_GVpair,%uni_GVpair,%varAAc_varint);
my $k=0;
for(my $i=0;$i<$#interactome+1;$i++){
  my @cmap = split(/\t/,$interactome[$i]);
#data in following format
#UniprotAC; pdb id; interface #; domain 1; PFAM id 1; chain id 1; shift 1; domain s 1; domain e 1; domain 2; PFAM id 2; chain id 2; shift 2; domain s 2; domain e 2; AA 1; AA id 1; AA 2; AA id 2; type; score; zscore

# get pdb and chain IDs and get uniprot annotation 
  my $pdbch1=$cmap[1].$cmap[5];
  my $pdbch2=$cmap[1].$cmap[11];
  my $shift1=$cmap[6];#difference between uniprot and pdb sequence numbering of partner 1
  my $shift2=$cmap[12];#difference between uniprot and pdb sequence numbering of partner 2
  my $uni=$cmap[0];#uniprotAC
  my $interface=$uni."_".$pdbch1."-".$pdbch2."-".$cmap[2];
  my $intfAA1=$uni."_".$pdbch1."_".$cmap[15];
  my $intfAA2=$uni."_".$pdbch2."_".$cmap[17];

#map variations and count occurences of contacts in the interface or in both faces separately
  for(my $i=0;$i<$#variat+1;$i++){
    my $waryjat=$variat[$i][6]."_".$variat[$i][2].$variat[$i][3].$variat[$i][4];

    if($variat[$i][6] eq $uni){#variant matches given uniprotID/entrezID
      my $pos = $variat[$i][3];
      my $intpos1 = $cmap[15]+$shift1;
      my $intpos2 = $cmap[17]+$shift2;
      if($pos == $intpos1) {# if variant is in contact (varAA number is same as AA1 in cmap of specific interface)
	if ($cmap[16] ne $variat[$i][3]){
	   print "WARNING: reference AA does not match interactome\n";
	   if ($cmap[16] ne $variat[$i][5]){print "WARNING: variant AA does not match interactome too\n";}
	   if ($cmap[16] eq $variat[$i][5]){print "WARNING: but variant AA match interactome\n";}   
	}
	$cmapAAisGV{$intfAA1}=1;
	$varAAc{$waryjat}++;#count contacts for given AA
	if($cmap[19] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
	elsif($cmap[19] eq "ms" or $cmap[19] eq "sm"){$varAAc_sm{$waryjat}++;}
	elsif($cmap[19] eq "mm"){$varAAc_mm{$waryjat}++;}
	$intfc{$interface}++;#count number of SNPs within specific interface in general
	$intfc_1{$interface}++;#or separately for 1st site
      }
      elsif($pos eq $intpos2){# if variant is in contact (varAA number is same as AA2 (mirror) in cmap of specific interface)
	if ($cmap[18] ne $variat[$i][3]){
	   print "WARNING: reference AA does not match interactome\n";
	   if ($cmap[18] ne $variat[$i][5]){print "WARNING: variant AA does not match interactome too\n";}
	   if ($cmap[18] eq $variat[$i][5]){print "WARNING: but variant AA match interactome\n";}   
	}
	$cmapAAisGV{$intfAA2}=1;
	$varAAc{$waryjat}++;# count contacts for given AA
	if($cmap[19] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
	elsif($cmap[19] eq "ms" or $cmap[19] eq "sm"){$varAAc_sm{$waryjat}++;}
	elsif($cmap[19] eq "mm"){$varAAc_mm{$waryjat}++;}
	$intfc{$interface}++;# count number of SNPs within specific interface in general
	$intfc_2{$interface}++;# or separately for mirror site
      }
    }
  }
}


#second round of comparison
my $name1="Edge_attributes.txt";
my $name2="IIN_graph.txt";
open(EDG, ">$name1");
open(NET, ">$name2");
for(my $i=0;$i<$#interactome+1;$i++){
  my @cmap = split(/\t/,$interactome[$i]);
#data in following format
#UniprotAC; pdb id; interface #; domain 1; PFAM id 1; chain id 1; shift 1; domain s 1; domain e 1; domain 2; PFAM id 2; chain id 2; shift 2; domain s 2; domain e 2; AA 1; AA id 1; AA 2; AA id 2; type; score; zscore

# get pdb and chain IDs and get uniprot annotation 
  my $pdbch1=$cmap[1].$cmap[5];
  my $pdbch2=$cmap[1].$cmap[11];
  my $shift1=$cmap[6];#difference between uniprot and pdb sequence numbering of partner 1
  my $shift2=$cmap[12];#difference between uniprot and pdb sequence numbering of partner 2
  my $uni=$cmap[0];#uniprotAC
  my $interface=$uni."_".$pdbch1."-".$pdbch2."-".$cmap[2];
  my $node1=$uni."_".$cmap[4]."_".$pdbch1."-".$cmap[2];	
  my $node2=$uni."_".$cmap[11]."_".$pdbch2."-".$cmap[2];	
  my $intfAA1=$uni."_".$pdbch1."_".$cmap[15];
  my $intfAA2=$uni."_".$pdbch2."_".$cmap[17];

#map variations
  for(my $i=0;$i<$#variat+1;$i++){
    my $waryjat=$variat[$i][6]."_".$variat[$i][2].$variat[$i][3].$variat[$i][4];
    if($variat[$i][6] eq $uni){#variant matches given uniprotID/entrezID
      my $pos = $variat[$i][3];
      my $intpos1 = $cmap[15]+$shift1;
      my $intpos2 = $cmap[17]+$shift2;
      printf NET "$node1\t$node2\t->\n";
      printf EDG "0\n";	
      if($pos == $intpos1) {#  if variant is in contact (varAA number is same as 1 AA in cmap of specific interface)
	if($cmapAAisGV{$intfAA2} > 0){
	  $varAAc_GVpair{$waryjat}++;#count number of contacts of given AA with residues that are also variants
	  $uni_GVpair{$uni}++;#for given uniprotAC count number of contacts with residues that are also variant
	  printf EDG "3\n";	
	}
	else{printf EDG "1\n";}	
	if($intfc_2{$interface} > 0){
	  $varAAc_varint{$waryjat}++;#count contacts of given AA with residue that contains GV
	}
      }
      if($pos == $intpos2) {#  if variant is in contact (varAA number is same as 1 AA in cmap of specific interface)
	if($cmapAAisGV{$intfAA1} > 0){
	  $varAAc_GVpair{$waryjat}++;#count number of contacts of given AA with residues that are also variants
	  $uni_GVpair{$uni}++;#for given uniprotAC count number of contacts with residues that are also variant
	   printf EDG "3\n";	
	}
	else{printf EDG "2\n";}	
	if($intfc_1{$interface} > 0){
	  $varAAc_varint{$waryjat}++;#count contacts of given AA with residue that contains GV
	}
      }
    }
  }
}


sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
		
