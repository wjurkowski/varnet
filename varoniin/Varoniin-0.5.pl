#!/usr/bin/perl -w
use strict;
use warnings;

#Two files are needed: 1) list of variations 2) 3did interactome file
if ($#ARGV != 1) {die "Program used with parameters [list of variations] [interactome]\n";}

#TODO: update description and implement all below listed features
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
my(@variat, %gvtype, %disease, %choroba, %allvarAA, %uniSNV);
my(%AAc,%AAc_ss,%AAc_sm,%AAc_ms,%AAc_mm);#AA counts
my(%varAAinItf,%varAApaired,%varAAlinked,%varAAc,%varAAc_ss,%varAAc_sm,%varAAc_ms,%varAAc_mm);#GV AA counts
my(%varAAc_GVpair,%varAAc_GVint,%varAAc_pair,%varAAc_int);#GV AA counts in context
for(my $i=0;$i<$#variations+1;$i++){
	my @tab = split(/\t/,$variations[$i]);#GVkey, Symbol, Entrez,  Uniprot AC, Uniprot ID, Ref AA, Position, Var AA, disease
	push (@{$variat[$i]}, @tab);
	my $waryjat=$variat[$i][3]."_".$variat[$i][5].$variat[$i][6];
	#initialize hashes for printing not mapped variants
	$varAAinItf{$waryjat}=0;
	$varAApaired{$waryjat}=0;
	$varAAlinked{$waryjat}=0;
	$varAAc{$waryjat}=0; 
	$varAAc_ss{$waryjat}=0;
	$varAAc_sm{$waryjat}=0;
	$varAAc_ms{$waryjat}=0;
	$varAAc_mm{$waryjat}=0;
	$varAAc_GVpair{$waryjat}=0;
	$varAAc_GVint{$waryjat}=0;
	$varAAc_pair{$waryjat}=0;
	$varAAc_int{$waryjat}=0;
	$gvtype{$waryjat}=$variat[$i][9];
	$variat[$i][8] =~ s/\s+$//;
	$variat[$i][8] =~ s/^\s+//;
	if(exists $disease{$waryjat}){#hash mapping all variants to disease
	    my $old= $disease{$waryjat};
	    my @test=split(";",$old);
	    unless (grep { $_ eq $variat[$i][8]} @test ){
	      my $new = $old.";".$variat[$i][8];
	      $disease{$waryjat}=$new;
	    }   
	}
	else{
	    $disease{$waryjat}=$variat[$i][8];
	}
	$allvarAA{$waryjat}=$i;
	$uniSNV{$variat[$i][3]}++;#hash counting number of GV in protein
#print "cece $waryjat\n";
#print "$variat[$i][3] $uniSNV{$variat[$i][3]}\n";
}	

#output files
open(NET, ">SNVmap-IIN_graph.txt");
open(NODE, ">SNVmap-IIN_nodes.txt");
open(IOUT, ">SNVmap-IIN.out");
open(POUT, ">SNVmap-PIN.out");
open(CTd, "SNVmap-CTdisease");
open(CTn, "SNVmap-CTnode");
open(RES, ">SNVmap-AAstats.out");
open(VARS, ">SNVmap-SNVstats.out");

#read in interactome: 3did, Jail? etc. 
#first round of input analysis
my @interactome=open_file($ARGV[1]);
my(%contacts,%residues);
my(%intfc,%ifc,%intfcgv,%ifcgv);#interface counts
my (%uniAAcint,%uniSNVcint,%edge);
my $k=0;

for(my $i=0;$i<$#interactome+1;$i++){
my @cmap = split(/\t/,$interactome[$i]);
#data in following format
#pdbid; interface#; uniprotAC1; domain1; PFAMid1; chid1; shift1; domain_s1; domain_e1; uniprotAC2; domain2; PFAMid2; chid2; shift2; domain_s2; domain_e2; AA1; AAid1; AA2; AAid2; type; score; zscore

# get pdb and chain IDs and get uniprot annotation 
  my $uni1=$cmap[2];#uniprotAC
  my $uni2=$cmap[9];#uniprotAC
  my $shift1=$cmap[6];#difference between uniprot and pdb sequence numbering of partner 1
  my $shift2=$cmap[13];#difference between uniprot and pdb sequence numbering of partner 2
  my $domain1=$cmap[3];#domain name
  my $domain2=$cmap[10];#domain name
  my $aar1=$cmap[7]."-".$cmap[8];#residue range
  my $aar2=$cmap[14]."-".$cmap[15];#residue range
  my $node1=$uni1."_".$domain1."_".$aar1;
  my $node2=$uni2."_".$domain2."_".$aar2;
  my $interface=$node1."<->".$node2;#definition of the interface
  my $resnAA1=$cmap[16]+$shift1;
  my $resnAA2=$cmap[18]+$shift2;
  my $residAA1=$cmap[17];
  my $residAA2=$cmap[19];
  my $reskey1=$uni1."_".$residAA1.$resnAA1;
  my $reskey2=$uni2."_".$residAA2.$resnAA2;
  my $key=$reskey1."-".$reskey2;
  
  $contacts{$key}=$interface;#link interacting residues with interface
  $residues{$reskey1}=$node1; 	
  $residues{$reskey2}=$node2; 	 	
#count all contacts and set variants to zero
#%AAc - total number of contacts of non SNV residue 
#%AAc_ss - for given residue number of contacts between sidechains
#%AAc_ms - for given residue number of contacts between its backbone and mirror sidechain
#%AAc_sm - for given residue number of contacts between its sidechain and mirror backbone
#%AAc_mm - for given residue number of contacts between backbones
#%varAAc - total number of contacts of given GV
#%varAAc_ss - for given GV number of contacts between sidechains
#%varAAc_sm - for given GV number of contacts between sidechain and backbone (or vice versa)
#%varAAc_mm - for given GV number of contacts between backbones
#%intfc - total number of contacts on interface
#%intfcgv - total number of SNV contacts on interface
#%ifc - number of contacts on faces/nodes
#%ifcgv - number of SNV contacts on faces/nodes
#%varAAc_GVpair - for given GV count contacts with other GV
#%varAAc_GVint - for given GV count contacts with interface carrying other GV
#%uniSNVcint - for given uniprotAC count contacts with interfaces carrying GV
#%varAAc_pair - for given variant count contacts with residues that are not variants 
#%varAAc_int - for given variant count contacts with interface NOT carrying other GV
#%uniAAcint - for given uniprotAC count contacts within interfaces not carrying GV

  $ifcgv{$node1}=0 unless (exists $ifcgv{$node1});#left
  $ifcgv{$node2}=0 unless (exists $ifcgv{$node2});#right
  $intfcgv{$interface}=0 unless (exists $intfcgv{$interface});
  #number of non SNV residues within interface
  #$intfc{$interface}=0;
  $uniAAcint{$uni1}=0 unless (exists $uniAAcint{$uni1});
  $uniAAcint{$uni2}=0 unless (exists $uniAAcint{$uni2});
  $uniSNVcint{$uni1}=0 unless (exists $uniSNVcint{$uni1});
  $uniSNVcint{$uni2}=0 unless (exists $uniSNVcint{$uni2});
  $disease{$reskey1}="NA" unless (exists $disease{$reskey1});
  $disease{$reskey2}="NA" unless (exists $disease{$reskey2}); 
  $gvtype{$reskey1}="NA" unless (exists $gvtype{$reskey1});
  $gvtype{$reskey2}="NA" unless (exists $gvtype{$reskey2});
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
  $choroba{$node1}="NA";
  $choroba{$node2}="NA";
  $intfc{$interface}++;
  $ifc{$node1}++;
  $ifc{$node2}++;

#check targets
if(exists $uniSNV{$uni1} or exists $uniSNV{$uni2}){
#print "$uni1 $uniSNV{$uni1} ground control interactome iteration number: $i\n";

#map variations and count occurences of contacts in the interface or in both faces separately
#for given uniprot AC of variant compare with uniprot AC of interactome	
#if uniprotAC and positions (after shift according to mapping) are the same map variant on interacting pair of residues 		
#count occurences of GV in interfaces etc.
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
		if($cmap[20] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
		elsif($cmap[20] eq "ms"){$varAAc_ms{$waryjat}++;}
		elsif($cmap[20] eq "sm"){$varAAc_sm{$waryjat}++;}
		elsif($cmap[20] eq "mm"){$varAAc_mm{$waryjat}++;}
		$intfcgv{$interface}++;#count number of SNPs within specific interface in general
		$ifcgv{$node1}++;#or separately for 1st site
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
		if($cmap[20] eq "ss"){$varAAc_ss{$waryjat}++;}#count contacts classified into groups (ss, ms, mm)
		elsif($cmap[20] eq "ms"){$varAAc_sm{$waryjat}++;}
		elsif($cmap[20] eq "sm"){$varAAc_ms{$waryjat}++;}
		elsif($cmap[20] eq "mm"){$varAAc_mm{$waryjat}++;}
		$intfcgv{$interface}++;# count number of SNPs within specific interface in general
		$ifcgv{$node2}++;# or separately for mirror site
	}
    }
  }
}
#else{print "REMARK protein not present on variants list\n";}
}

#second round of comparison
#analyse all contacting residues identified by combination of reside identifiers
#counts contacts of SNV in mirror face and interface context
foreach my $key (keys %contacts){#iterate all contacts hashed by combination of residue identifiers (uniprot-AAposAAid)
	my @klucze=split("-",$key);
	my $reskey1=$klucze[0];
	my $reskey2=$klucze[1];
	my @tem1=split("_",$reskey1);
	my $uni1=$tem1[0];
	@tem1=split("_",$reskey2);
	my $uni2=$tem1[0];
	my $intf=$contacts{$key};
	# get interface and split into graph nodes
	my @faces=split("<->",$intf);
	my $node1=$faces[0];
	my $node2=$faces[1];
	$uniSNV{$uni1}=0 unless (exists $uniSNV{$uni1});
	$uniSNV{$uni2}=0 unless (exists $uniSNV{$uni2});

	if(exists $varAAc{$reskey1}){
		$choroba{$node1}=$disease{$reskey1};
		$varAAinItf{$reskey1}=1;
		if($varAAc{$reskey2}){# mirror residue is also GV
			$varAApaired{$reskey1}=1;
			$varAApaired{$reskey2}=1;
			$varAAc_GVpair{$reskey1}++;
			$varAAc_GVpair{$reskey2}++;
			$edge{$intf}=3;
			$varAAc_GVint{$reskey1}++;
			$varAAc_GVint{$reskey2}++;	
			$uniSNVcint{$uni1}++;
			$uniSNVcint{$uni2}++;	
		}
		else{# mirror residue is not GV
			$varAAc_pair{$reskey1}++;
			$edge{$intf}=1 unless($edge{$intf} > 0);
			if($intfc{$intf} >1){# but interface contains other GV
				$varAAlinked{$reskey1}=1;
				$varAAc_GVint{$reskey1}++;
				$uniSNVcint{$uni1}++;
				$uniSNVcint{$uni2}++;
				$edge{$intf}=2 unless($edge{$intf} > 1);
			}
			else{# interface does not contain other GV 
				$varAAc_int{$reskey1}++;
				$uniAAcint{$uni1}++;
				$uniAAcint{$uni2}++;
			}
		}
	}
	elsif(exists $varAAc{$reskey2}){
		$choroba{$node2}=$disease{$reskey2};
		$varAAinItf{$reskey2}=1;
		unless(exists $varAAc{$reskey1}){# mirror residue is not GV
			$varAAc_pair{$reskey2}++;
			$edge{$intf}=1 unless($edge{$intf} > 0);
			if($intfc{$intf} >1){# but interface contains other GV
				$varAAlinked{$reskey2}=1;
				$varAAc_GVint{$reskey2}++;
				$uniSNVcint{$uni1}++;
				$uniSNVcint{$uni2}++;
				$edge{$intf}=2 unless($edge{$intf} > 1);
			}
			else{# interface does not contain other GV
				$varAAc_int{$reskey2}++;
				$uniAAcint{$uni1}++;
				$uniAAcint{$uni2}++;
			}
		}
	}
}

my (%innodeSNVload,%innodeSNVdisease,%outnodeSNVload,%outnodeSNVdisease,$SNVin,$SNVout);
#contact map stats: number of contacts of different types calculated for single amino acid classified as reference or GV
#analyze and print per residue'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''..
printf "variantAA\tdisease\tgvtype\tAAc\tAAc_ss\tAAc_sm\tAAc_ms\tAAc_mm\n";
foreach my $reskey (keys %AAc){
	printf RES "$reskey\t$disease{$reskey}\t$gvtype{$reskey}\t$AAc{$reskey}\t$AAc_ss{$reskey}\t$AAc_sm{$reskey}\t$AAc_ms{$reskey}\t$AAc_mm{$reskey}\n";
	#define hashes storing number of SNV in and outside interface
	$innodeSNVload{$residues{$reskey}}=0;
	$innodeSNVdisease{$disease{$reskey}}=0;
	$outnodeSNVload{$residues{$reskey}}=0;
	$outnodeSNVdisease{$disease{$reskey}}=0;
}

#analyze and print per variant
#count SNV in and out of interfaces 
printf VARS "variantAA\tdisease\tgvtype\tvarAAinItf\tvarAApaired\tvarAAlinked\tvarAAc\tvarAAc_ss\tvarAAc_sm\tvarAAc_ms\tvarAAc_mm\tvarAAc_GVpair\tvarAAc_GVint\tvarAAc_pair\tvarAAc_int\n";
foreach my $reskey (keys %allvarAA){
	if($varAAinItf{$reskey} == 1){
		$innodeSNVload{$residues{$reskey}}++;#count SNV on a node
		$innodeSNVdisease{$disease{$reskey}}++;
		$SNVin++;	
	}
	elsif($varAAinItf{$reskey} == 0){
		$outnodeSNVload{$residues{$reskey}}++;#count SNV outside any node
		$outnodeSNVdisease{$disease{$reskey}}++;#count SNV outside any node for particular disease
		$SNVout++;
	}
	printf VARS "$reskey\t$disease{$reskey}\t$gvtype{$reskey}\t$varAAinItf{$reskey}\t$varAApaired{$reskey}\t$varAAlinked{$reskey}\t$varAAc{$reskey}\t$varAAc_ss{$reskey}\t$varAAc_sm{$reskey}\t$varAAc_ms{$reskey}\t$varAAc_mm{$reskey}\t$varAAc_GVpair{$reskey}\t$varAAc_GVint{$reskey}\t$varAAc_pair{$reskey}\t$varAAc_int{$reskey}\n";    
}

#Per disease contingency table of variants in and out of interfaces
printf CTd	"disease\tSNVinInterface\tSNVoutInterface\n";
foreach my $dis (keys %innodeSNVdisease){
	printf CTd "$dis\t$innodeSNVdisease{$dis}\t$outnodeSNVdisease{$dis}\n";
}

#Per node contingency table of variants in and out of interfaces
printf CTn "node\tSNVinInterface\tSNVoutInterface\n";
foreach my $node (keys %innodeSNVload){
	printf CTn "$node\t$innodeSNVload{$node}\t$outnodeSNVload{$node}\n";
}

print "Contingency table of total number of variants in and out of interfaces (compared for same nodes)";
print "SNVinInterface\tSNVoutInterface\n";
print "$SNVin\t$SNVout\n";

#analyze and print per node (half of interface)
#print graph and attributes
#print 
#edge: 0 - no variants, 1 - GV in contact, 2 - GV in interface carrying other GV, 3 - GV in contact with other GV
printf IOUT "Interface\tAAc\tSNVc\tdiseaseN1\tAAcNode1\tSNVcNode1\tdiseaseN2\tAAcNode2\tSNVcNode2\n";
foreach my $ifa (keys %intfc){
# get interface and split into faces
my @faces=split("<->",$ifa);
my $node1=$faces[0];
my @tem=split("_",$node1);
my $uni1 = $tem[0];
my $node2=$faces[1];
@tem=split("_",$node2);
my $uni2 = $tem[0];
	printf NET "$node1\t$node2\t->\t$node1 (pp) $node2\t$edge{$ifa}\t$intfcgv{$ifa}\n";	
	printf NODE "$node1\t$ifc{$node1}\t$choroba{$node1}\n";
	printf NODE "$node2\t$ifc{$node2}\t$choroba{$node2}\n";
#numbers of contacts in/out of interface 
	printf IOUT "$ifa\t$intfc{$ifa}\t$intfcgv{$ifa}\t$choroba{$node1}\t$ifc{$node1}\t$ifcgv{$node1}\t$choroba{$node2}\t$ifc{$node2}\t$ifcgv{$node2}\n";
}

#analyze and print per protein
# print number of GV in protein
foreach my $prot (keys %uniSNV){
	printf POUT "SNVinProt1\tSNVcontInProt1\tAAcontInProt1\tSNVinProt2\tSNVcontInProt2\tAAcontInProt2\n";
	printf POUT "$uniSNV{$prot}\t$uniSNVcint{$prot}\t$uniAAcint{$prot}\t$uniSNV{$prot}\t$uniSNVcint{$prot}\t$uniAAcint{$prot}\n";
}

sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
		
