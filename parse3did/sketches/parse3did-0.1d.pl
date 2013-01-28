#!/usr/bin/perl -w
use strict;
use warnings;

#Files needed on input: 3did interactome from 3did_flat
if ($#ARGV != 1) {die "Program used with parameters [interactome] [sifts uniprot-pdb mapping]\n";}


#parse PDB swissport mapping
my @swifts=open_file($ARGV[1]);
my (%pdb2uni);
for (my $i=0;$i<=$#swifts;$i++){
	my @lin=split(/\t/,$swifts[$i]);
	my $pdbch=$lin[0].$lin[1];;
	$pdb2uni{$pdbch}=$i;	
}

#get interfaces from 3did
my @threedid=open_file($ARGV[0]);
my $interactome="3did_interactome.tsv";
my $ddi="3did-DDI-graph.tsv";
my $iii="3did-III-graph.tsv";
open(OUT, ">$interactome");
open(OUT2, ">$ddi");
open(OUT3, ">$iii");
my $str="NO";
my $oldpdb="NULL";
my (@in,$uniprotAC1,$uniprotAC2,$shift1,$shift2,$ni,$face1,$face2);
#parse interactions from 
#Other releases will include other sources of interactions
printf OUT "pdbID\tintf#\tUP_AC1\tdomain 1\tPFAM id 1\tchID1\tshift1\td_s1\td_e1\tUP_AC2\tdomain 2\tPFAM id 2\tchID2\tshift2\td_s2\td_e2\tAA1\tAAid1\tAA2\tAAid2\ttype\tscore\tzscore\n";
printf OUT2 "UniprotAC1-domain1\tUniprotAC2-domain2\n";
printf OUT3 "UniprotAC1-domain-face\tUniprotAC2-domain-face\n";
for(my $i=0;$i<$#threedid;$i++){
	if($threedid[$i] eq "//"){next;}# skip line
	my @int = split(/\t/,$threedid[$i]);
	my $id = $int[0];
	$id =~ s/\W=//;	
	if($id eq "ID"){
		$in[0]=$int[1];#name domain 1
		$in[1]=$int[2];#name domain 2
		$in[2]=substr($int[3],0,index($int[3],"@"));#PFAM 1 id	
		$in[3]=substr($int[4],0,index($int[4],"@"));#PFAM 2 id
		$in[2] =~ s/\s*\(//;
		$in[3] =~ s/\s*\)//;
		$str="NO";
		$ni=1;
	} 
	elsif($id eq "3D"){
		$in[4]=$int[1];#pdb id
		if($in[4] eq $oldpdb){$ni++;}
		else{$ni=1;}	
		my @tab1=split(":",$int[2]);
		$in[5]=$tab1[0];#domain 1 chain id
		my @tab2=split("-",$tab1[1]);
		$in[6]=$tab2[0];#domain 1 start 
		$in[7]=$tab2[1];#domain 1 end
		my @tab3=split(":",$int[3]);
		$in[8]=$tab3[0];#domain 2 chain id
		my @tab4=split("-",$tab3[1]);
		$in[9]=$tab4[0];#domain 2 start
		$in[10]=$tab4[1];#domain 2 end
		$in[11]=$int[4];#score
		$in[12]=$int[5];#zscore
		$str="YES";		
		my $pdbch1=$in[4].$in[5];
		my $pdbch2=$in[4].$in[8];
		if(exists $pdb2uni{$pdbch1}){	
			my @mapping=split(/\t/,$swifts[$pdb2uni{$pdbch1}]);
			$uniprotAC1=$mapping[2];
			my $UPseq_s=$mapping[7];
			my $PDBseq_s=$mapping[5];
			if($PDBseq_s =~ /[A-Z]/){$PDBseq_s=clean_resn($PDBseq_s);}
			$shift1=$UPseq_s-$PDBseq_s;
		}
		else{
			$uniprotAC1="None";
			$shift1=0;
		}
		if(exists $pdb2uni{$pdbch2}){	
		my @mapping=split(/\t/,$swifts[$pdb2uni{$pdbch2}]);
			$uniprotAC2=$mapping[2];
			my $UPseq_s=$mapping[7];
			my $PDBseq_s=$mapping[5];
			if($PDBseq_s =~ /[A-Z]/){$PDBseq_s=clean_resn($PDBseq_s);}
			$shift2=$UPseq_s-$PDBseq_s;
		}
		else{
			$uniprotAC2="None";
			$shift2=0;
		}
		$oldpdb=$in[4];
		
		if($in[6] =~ /[A-Z]/){$in[6]=clean_resn($in[6]);}
		if($in[7] =~ /[A-Z]/){$in[7]=clean_resn($in[7]);}
		if($in[9] =~ /[A-Z]/){$in[9]=clean_resn($in[9]);}
		if($in[10] =~ /[A-Z]/){$in[10]=clean_resn($in[10]);}
		my $d1s=$in[6]+$shift1;
print "$int[2]\n";			
		my $d1e=$in[7]+$shift1;
		my $d2s=$in[9]+$shift2;
		my $d2e=$in[10]+$shift2;
		$face1=$d1s."-".$d1e;
		$face2=$d2s."-".$d2e;
	}	
	else{
	 if($str eq "YES"){
		my $node1=$uniprotAC1."-".$in[0]."-".$face1;
		my $node2=$uniprotAC2."-".$in[1]."-".$face2;
		#print interacting AA pairs in given interface. 
		#Interface is defined as instance of pfam domain pairs in contact in particular PDB file. Given PDB file may include multiple interfaces representing same or different pfam-pfam interaction
		#Interface # is to identify same pfam-pfam in unique pdb
		#pdb id; interface #; $uniprotAC1 domain 1; PFAM id 1; chain id 1; shift 1; domain s 1; domain e 1; $uniprotAC2; domain 2; PFAM id 2; chain id 2; shift 2; domain s 2; domain e 2; AA 1; AA id 1; AA 2; AA id 2; type; score; zscore
		printf OUT "$in[4]\t$ni\t$uniprotAC1\t$in[0]\t$in[2]\t$in[5]\t$shift1\t$in[6]\t$in[7]\t$uniprotAC2\t$in[1]\t$in[3]\t$in[8]\t$shift2\t$in[9]\t$in[10]\t$int[2]\t$int[0]\t$int[3]\t$int[1]\t$int[4]\t$in[11]\t$in[12]\n";
		printf OUT2 "$uniprotAC1-$in[0]\t$uniprotAC2-$in[1]\n";
		printf OUT3 "$node1\t$node2\n";
	 }
	}
}

sub clean_resn{
	my ($dirty)=@_;
	my $letter=$dirty;
	$letter=~s/\d//g;
	my $value=$dirty;
	$value=~s/\D//g;
	my $numb=ord($letter);
    my $clean=$value-$numb;
	return $clean;
}

sub open_file{
        my ($file_name)=@_;
        open(INP1, "< $file_name") or die "Can not open an input file: $!";
        my @file1=<INP1>;
        close (INP1);
        chomp @file1;
        return @file1;
}		
