#!/usr/bin/perl -w
use strict;
use warnings;

#Two files are needed: 1) list of variations 2) 3did interactome file
if ($#ARGV != 1) {die "Program used with parameters [list of variations] [interactome]\n";}

#TODO: update description and implement all elow listed features
#query interface from 3did
#	print structuraly resolved interactions
#	interpret interfaces
#		pfam, functional residues, complexation residues		
#read pfam - pdb mapping from pfam website: pfam definitions according to pdb resid numbering

#from list of variations
	#map variations on the interactome graph
	#identify pairs or chains of interactions = > group 
#from list of proteins
#	get interactions
#	get domain pairs (of interacting protein pair)
#	get list of interacting AA

#parse_variations
my @variations=open_file($ARGV[0]);
my(@variat, %allvar, %unic);
my(%AAc,%AAc_ss,%AAc_sm,%AAc_ms,%AAc_mm,%varAAc,%varAAout,%varAAc_ss,%varAAc_sm,%varAAc_ms,%varAAc_mm);#AA counts
my (%varAAc_GVpair,%varAAc_GVint,%varAAc_pair,%varAAc_int);#AA counts in context
for(my $i=0;$i<$#variations+1;$i++){
	my @tab = split(/\t/,$variations[$i]);#GVkey, Symbol, Entrez,  Uniprot AC, Uniprot ID, Ref AA, Position, Var AA, disease
	push (@{$variat[$i]}, @tab);
	my $waryjat=$variat[$i][3]."_".$variat[$i][5].$variat[$i][6];
	#initialize hashes for printing not mapped variants
	$varAAc{$waryjat}=0; 
	$varAAout{$waryjat}=0;
	$varAAc_ss{$waryjat}=0;
	$varAAc_sm{$waryjat}=0;
	$varAAc_ms{$waryjat}=0;
	$varAAc_mm{$waryjat}=0;
	$varAAc_GVpair{$waryjat}=0;
	$varAAc_GVint{$waryjat}=0;
	$varAAc_pair{$waryjat}=0;
	$varAAc_int{$waryjat}=0;
	$allvar{$waryjat}=$variat[$i][8];#hash mapping all variants to disease
	$unic{$variat[$i][3]}++;#hash counting number of GV in protein
#print "cece $waryjat\n";
#print "$variat[$i][3] $unic{$variat[$i][3]}\n";
}	

#output files
my $name="SNVmapped-IIN_graph.txt";
my $name2="SNVmapped-IIN_nodes.txt";
my $name3="SNVmapped-IIN.out";
my $name4="AAstats.out";
my $name5="GVstats.out";
open(NET, ">$name");
open(NODE, ">$name2");
open(GOUT, ">$name3");
open(RES, ">$name4");
open(VARS, ">$name5");

#read in interactome
my @interactome=open_file($ARGV[1]);
my(%contacts,%disease);
my(%intfc,%ifc,%intfgvc,%intfgvc1,%intfgvc2);#interface counts
my (%uni_int,%uni_GVint,%edge);
my $k=0;

for(my $i=0;$i<$#interactome+1;$i++){
  my @cmap = split(/\t/,$interactome[$i]);
#data in following format
#pdbid; interface#; uniprotAC1; domain1; PFAMid1; chid1; shift1; domain_s1; domain_e1; uniprotAC2; domain2; PFAMid2; chid2; shift2; domain_s2; domain_e2; AA1; AAid1; AA2; AAid2; type; score; zscore

# get pdb and chain IDs and get uniprot annotation 
  my $shift1=$cmap[6];#difference between uniprot and pdb sequence numbering of partner 1
  my $shift2=$cmap[13];#difference between uniprot and pdb sequence numbering of partner 2
  my $uni1=$cmap[2];#uniprotAC
  my $uni2=$cmap[9];#uniprotAC
  my $domain1=$cmap[3];#domain name
  my $domain2=$cmap[10];#domain name
  my $aar1=$cmap[7]."-".$cmap[8];#residue range
  my $aar2=$cmap[14]."-".$cmap[15];#residue range
  my $face1=$uni1."_".$domain1."_".$aar1;
  my $face2=$uni2."_".$domain2."_".$aar2;
  my $interface=$face1."<->".$face2;#definition of the interface
  my $resnAA1=$cmap[16]+$shift1;
  my $resnAA2=$cmap[18]+$shift2;
  my $residAA1=$cmap[17];
  my $residAA2=$cmap[19];
  my $reskey1=$uni1."_".$residAA1.$resnAA1;
  my $reskey2=$uni2."_".$residAA2.$resnAA2;
  my $key=$reskey1."-".$reskey2;
  
  $intfgvc1{$interface}=0 unless (exists $intfgvc1{$interface});#left
  $intfgvc2{$interface}=0 unless (exists $intfgvc2{$interface});#right
  $intfgvc{$interface}=0 unless (exists $intfgvc{$interface});
  #number of non SNV residues within interface
  #$intfc{$interface}=0;
  
  $contacts{$key}=$interface;#link interacting residues with interface

#count all contacts and set variants to zero
#%AAc - total number of contacts of non SNV residue 
#%AAc_ss - for given residue number of contacts between sidechains
#%AAc_ms - for given residue number of contacts between its backbone and mirror sidechain
#%AAc_sm - for given residue number of contacts between its sidechain and mirror backbone
#%AAc_mm - for given residue number of contacts between backbones
#%intfc - total number of non GV in interface
#%intfc1 - number of non GV in left face
#%intfc2 - number of GV in right face
  $unic{$uni1}=0 unless (exists $unic{$uni1});
  $unic{$uni2}=0 unless (exists $unic{$uni2});
  $uni_int{$uni1}=0 unless (exists $uni_int{$uni1});
  $uni_int{$uni2}=0 unless (exists $uni_int{$uni2});
  $uni_GVint{$uni1}=0 unless (exists $uni_GVint{$uni1});
  $uni_GVint{$uni2}=0 unless (exists $uni_GVint{$uni2});
  $disease{$reskey1}="NA" unless (exists $disease{$reskey1});
  $disease{$reskey2}="NA" unless (exists $disease{$reskey2});
  $varAAc{$reskey1}=0 unless (exists $varAAc{$reskey1});
  $varAAc{$reskey2}=0 unless (exists $varAAc{$reskey2});
  $varAAc_ss{$reskey1}=0 unless (exists $varAAc_ss{$reskey1});
  $varAAc_ss{$reskey2}=0 unless (exists $varAAc_ss{$reskey2});
  $varAAc_ms{$reskey1}=0 unless (exists $varAAc_ms{$reskey1});
  $varAAc_ms{$reskey2}=0 unless (exists $varAAc_ms{$reskey2});
  $varAAc_sm{$reskey1}=0 unless (exists $varAAc_sm{$reskey1});
  $varAAc_sm{$reskey2}=0 unless (exists $varAAc_sm{$reskey2});
  $varAAc_mm{$reskey1}=0 unless (exists $varAAc_mm{$reskey1});
  $varAAc_mm{$reskey2}=0 unless (exists $varAAc_mm{$reskey2});
  $AAc_ss{$reskey1}=0 unless (exists $AAc_ss{$reskey1});
  $AAc_ss{$reskey2}=0 unless (exists $AAc_ss{$reskey2});
  $AAc_ms{$reskey1}=0 unless (exists $AAc_ms{$reskey1});
  $AAc_ms{$reskey2}=0 unless (exists $AAc_ms{$reskey2});
  $AAc_sm{$reskey1}=0 unless (exists $AAc_sm{$reskey1});
  $AAc_sm{$reskey2}=0 unless (exists $AAc_sm{$reskey2});
  $AAc_mm{$reskey1}=0 unless (exists $AAc_mm{$reskey1});
  $AAc_mm{$reskey2}=0 unless (exists $AAc_mm{$reskey2});
  $AAc{$reskey1}++;# count contacts for given AA
  $AAc{$reskey2}++;# count contacts for given AA
  if($cmap[20] eq "ss"){#count contacts classified into groups (ss, ms, mm)
    $AAc_ss{$reskey1}++;
    $AAc_ss{$reskey2}++;
  }
  elsif($cmap[20] eq "ms"){
    $AAc_ms{$reskey1}++;
    $AAc_sm{$reskey2}++;
  }
  elsif($cmap[20] eq "sm"){
    $AAc_ms{$reskey2}++;
    $AAc_sm{$reskey1}++;
  }
  elsif($cmap[20] eq "mm"){
    $AAc_mm{$reskey1}++;
    $AAc_mm{$reskey2}++;
  }
  $edge{$interface}=0 unless (exists  $edge{$interface});
  $intfc{$interface}++;
  $ifc{$face1}++;
  $ifc{$face2}++;
#map variations and count occurences of contacts in the interface or in both faces separately
#for given uniprot AC of variant compare with uniprot AC of interactome	
#if uniprotAC and positions (after shift according to mapping) are the same map variant on interacting pair of residues 		
#count occurences of GV in interfaces etc.
#
#%varAAgvc - total number of contacts of given GV
#%varAAgvc_ss - for given GV number of contacts between sidechains
#%varAAgvc_sm - for given GV number of contacts between sidechain and backbone (or vice versa)
#%varAAgvc_mm - for given GV number of contacts between backbones
#%intfgvc - total number of GV in interface
#%intfgvc1 - number of GV in left face
#%intfgvc2 - number of GV in right face
  for(my $i=0;$i<$#variat+1;$i++){
    my $waryjat=$variat[$i][3]."_".$variat[$i][5].$variat[$i][6];
    if($variat[$i][3] eq $uni1){#variant matches uniprotAC of face1
	my $pos = $variat[$i][6];
	if($pos == $resnAA1) {# if variant is in contact (varAA number is same as AA1 in cmap of specific interface)
		if ($residAA1 ne $variat[$i][5]){
		  print "WARNING: reference AA does not match interactome\n";
		  if ($residAA1 ne $variat[$i][7]){print "WARNING: variant AA does not match interactome too\n";}
		  if ($residAA1 eq $variat[$i][7]){print "WARNING: but variant AA match interactome\n";}   
		}
		$varAAc{$waryjat}++;#count contacts for given GV
		$disease{$waryjat}=$variat[$i][8];
		if($cmap[20] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
		elsif($cmap[20] eq "ms"){$varAAc_ms{$waryjat}++;}
		elsif($cmap[20] eq "sm"){$varAAc_sm{$waryjat}++;}
		elsif($cmap[20] eq "mm"){$varAAc_mm{$waryjat}++;}
		$intfgvc{$interface}++;#count number of SNPs within specific interface in general
		$intfgvc1{$interface}++;#or separately for 1st site
	}
	else{
		$varAAout{$waryjat}++;
	}
    }	
    elsif($variat[$i][3] eq $uni2){#variant matches uniprotAC of face2	
	my $pos = $variat[$i][6];
	if($pos eq $resnAA2){# if variant is in contact (varAA number is same as AA2 (mirror) in cmap of specific interface)
		if ($residAA2 ne $variat[$i][7]){
			print "WARNING: reference AA does not match interactome\n";
			if ($residAA2 ne $variat[$i][5]){print "WARNING: variant AA does not match interactome too\n";}
			if ($residAA2 eq $variat[$i][5]){print "WARNING: but variant AA match interactome\n";}   
		}	
		$varAAc{$waryjat}++;# count contacts for given AA
		$disease{$waryjat}=$variat[$i][8];
		if($cmap[20] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
		elsif($cmap[20] eq "ms"){$varAAc_sm{$waryjat}++;}
		elsif($cmap[20] eq "sm"){$varAAc_ms{$waryjat}++;}
		elsif($cmap[20] eq "mm"){$varAAc_mm{$waryjat}++;}
		$intfgvc{$interface}++;# count number of SNPs within specific interface in general
		$intfgvc2{$interface}++;# or separately for mirror site
	}
	else{
		$varAAout{$waryjat}++;
	}
    }
  }
}

#second round of comparison
#analyse all contacting residues identified by combination of reside identifiers
#%varAAc_GVpair - for given GV count contacts with other GV
#%varAAc_GVint - for given GV count contacts with interface carrying other GV
#%uni_GVint - for given uniprotAC count contacts with interfaces carrying GV
#%varAAc_pair - for given variant count contacts with residues that are not variants 
#%varAAc_int - for given variant count contacts with interface NOT carrying other GV
#%uni_int - for given uniprotAC count contacts within interfaces not carrying GV

#get GV counts in mirror face and interface context
foreach my $key (keys %contacts){#iterate all contacts hashed by combination of residue identifiers (uniprot-AAposAAid)
	my @klucze=split("-",$key);
	my $pol1=$klucze[0];
	my $pol2=$klucze[1];
	my @tem1=split("_",$pol1);
	my $uni1=$tem1[0];
	@tem1=split("_",$pol2);
	my $uni2=$tem1[0];
	my $intf=$contacts{$key};
	# get interface and split into graph nodes
	my @faces=split("<->",$intf);
	my $node1=$faces[0];
	my $node2=$faces[1];

	if(exists $varAAc{$pol1}){
		if($varAAc{$pol2}){# mirror residue is also GV
			$varAAc_GVpair{$pol1}++;
			$varAAc_GVpair{$pol2}++;
			$edge{$intf}=3;
			$varAAc_GVint{$pol1}++;
			$varAAc_GVint{$pol2}++;	
			$uni_GVint{$uni1}++;
			$uni_GVint{$uni2}++;	
		}
		else{# mirror residue is not GV
			$varAAc_pair{$pol1}++;
			$edge{$intf}=1 unless($edge{$intf} > 0);
			if($intfc{$intf} >1){# but interface contains other GV
				$varAAc_GVint{$pol1}++;
				$uni_GVint{$uni1}++;
				$uni_GVint{$uni2}++;
				$edge{$intf}=2 unless($edge{$intf} > 1);
			}
			else{# interface does not contain other GV 
				$varAAc_int{$pol1}++;
				$uni_int{$uni1}++;
				$uni_int{$uni2}++;
			}
		}
	}
	elsif(exists $varAAc{$pol2}){
		unless(exists $varAAc{$pol1}){# mirror residue is not GV
			$varAAc_pair{$pol2}++;
			$edge{$intf}=1 unless($edge{$intf} > 0);
			if($intfc{$intf} >1){# but interface contains other GV
				$varAAc_GVint{$pol2}++;
				$uni_GVint{$uni1}++;
				$uni_GVint{$uni2}++;
				$edge{$intf}=2 unless($edge{$intf} > 1);
			}
			else{# interface does not contain other GV
				$varAAc_int{$pol2}++;
				$uni_int{$uni1}++;
				$uni_int{$uni2}++;
			}
		}
	}
}

#analyze and print per residue
foreach my $face (keys %AAc){
	printf RES "$face\t$disease{$face}\t$AAc{$face}\t$varAAc{$face}\t$AAc_ss{$face}\t$varAAc_ss{$face}\t$AAc_sm{$face}\t$varAAc_sm{$face}\t$AAc_ms{$face}\t$varAAc_ms{$face}\t$AAc_mm{$face}\t$varAAc_mm{$face}\n";
}
#analyze and print per variant
foreach my $face (keys %allvar){
	printf VARS "$face\t$allvar{$face}\t$varAAc{$face}\t$varAAc_ss{$face}\t$varAAc_sm{$face}\t$varAAc_ms{$face}\t$varAAc_mm{$face}\t$varAAc_GVpair{$face}\t$varAAc_GVint{$face}\t$varAAc_pair{$face}\t$varAAc_int{$face}\n";    
}

#analyze and print per protein
# print number of GV in protein
#print graph and attributes
#print 
#edge: 0 - no variants, 1 - GV in contact, 2 - GV in interface carrying other GV, 3 - GV in contact with other GV
foreach my $ifa (keys %intfc){
# get interface and split into faces
my @faces=split("<->",$ifa);
my $node1=$faces[0];
my @tem=split("_",$node1);
my $uni1 = $tem[0];
my $node2=$faces[1];
@tem=split("_",$node2);
my $uni2 = $tem[0];
	printf NET "$node1\t$node2\t->\t$node1 (pp) $node2\t$edge{$ifa}\t$intfgvc{$ifa}\n";	
	printf NODE "$node1\t$ifc{$node1}\n";
	printf NODE "$node2\t$ifc{$node2}\n";
	printf GOUT "$ifa\t$intfc{$ifa}\t$intfgvc{$ifa}\t$ifc{$node1}\t$intfgvc1{$ifa}\t$ifc{$node2}\t$intfgvc2{$ifa}\t$unic{$uni1}\t$uni_GVint{$uni1}\t$uni_int{$uni1}\t$unic{$uni2}\t$uni_GVint{$uni2}\t$uni_int{$uni2}\n";
}


sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
		
