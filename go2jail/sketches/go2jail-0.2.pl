#!/usr/bin/perl -w
use strict;
use warnings;

#Files needed on input: jail interactoms
#if ($#ARGV != 2) {die "Program used with parameters [interactome] [sifts uniprot-pdb mapping] [pfam mapping]\n";}
if ($#ARGV != 1) {die "Program used with parameters [jail] [sifts]\
jail - directory with jail README files (http://bioinf-services.charite.de/jail/)
sifts - uniprotAC-pdbID mapping file (http://www.ebi.ac.uk/pdbe/docs/sifts/)
\n";}

#parse PDB swissport mapping
my @swifts=open_file($ARGV[1]);
my (%pdb2uni);
for (my $i=0;$i<=$#swifts;$i++){
	my @lin=split(/\t/,$swifts[$i]);
	my $pdbch=$lin[0].$lin[1];;
	$pdb2uni{$pdbch}=$i;	
}

#parse PDB Pfam mapping
#my @pfams=open_file($ARGV[1]);
#my (%pdb2pfam);
#for (my $i=0;$i<=$#pfams;$i++){
#	my @lin=split(/\t/,$pfams[$i]);
#	my $pdbch=$lin[0].$lin[1];;
#	$pdb2uni{$pdbch}=$i;	
#2PFL    A       22      613     PF02901.10      PFL
#}

#get interfaces from jail
my $jin1="$ARGV[0]/README_BIOUNIT_95_FAMILY_REDUCED.txt";
my $jin2="$ARGV[0]/README_CHAINS_95_FAMILY_REDUCED.txt";
my $jin3="$ARGV[0]/README_SCOP_FAMS.txt";
my $jin4="$ARGV[0]/README_NUCLEICS_95.txt";
my @jail1=open_file($jin1);
my @jail2=open_file($jin2);
my @jail3=open_file($jin3);
my @jail4=open_file($jin4);
my $interactome="jail_interactome.tsv";
my $ddi="jail-DDI-graph.tsv";
my $nucl="jail-na_binding.tsv";

open(OUT, ">$interactome");
open(OUT2, ">$ddi");
open(OUT3, ">$nucl");
#parse interactions from 
#Other releases will include other sources of interactions
printf OUT "pdbID\tchid1\tUP_AC1\tdomain1\tshift1\tchid2\tUP_AC2\tdomain2\tshift2\n";
printf OUT2 "UniprotAC1-domain1\tUniprotAC2-domain2\n";
printf OUT3 "pdbID\tchid1\tUP_AC1\tdomain1\tshift1\n";

get_jail1(@jail1);
get_jail1(@jail2);
get_jail2(@jail3);
get_jail3(@jail4);

sub get_jail1{
my (@jail)=@_;
for(my $i=0;$i<$#jail;$i++){
	my ($uniprotAC1,$uniprotAC2,$shift1,$shift2);
	if($jail[$i] =~ /^#/){next;}# skip line
	if($jail[$i] =~ /^$/){next;}# skip line
	my @int = split(/\t/,$jail[$i]);
	my $pdb	= $int[1];
	my $chid1 = $int[2];
	my $chid2 = $int[7];
	my $pdbch1=$pdb.$chid1;
print "cece $int[1]\n";
	my $pdbch2=$pdb.$chid2;
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
	printf OUT "$pdb\t$chid1\t$uniprotAC1\t$shift1\t$chid2\t$uniprotAC2\t$shift2\n";
	printf OUT2 "$uniprotAC1-$pdb-$chid1\t$uniprotAC2-$pdb-$chid2\n";
}
}

sub get_jail2{
my (@jail)=@_;
for(my $i=0;$i<$#jail;$i++){
	my ($uniprotAC1,$uniprotAC2,$shift1,$shift2);
	if($jail[$i] =~ /^#/){next;}# skip line
	if($jail[$i] =~ /^$/){next;}# skip line
	my @int = split(/\t/,$jail[$i]);
	my $scop1 = $int[1];
	my $scop2 = $int[8];
	my $pdb1= substr($scop1,1,4);
	my $pdb2= substr($scop2,1,4);
	my $chid1 = ucfirst(substr($scop1,5,1));
	my $chid2 = ucfirst(substr($scop2,5,1));
	my $pdbch1=$pdb1.$chid1;
	my $pdbch2=$pdb2.$chid2;
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
	printf OUT "$pdb1\t$chid1\t$uniprotAC1\t$shift1\t$pdb2\t$chid2\t$uniprotAC2\t$shift2\n";
	printf OUT2 "$uniprotAC1-$pdb1-$chid1\t$uniprotAC2-$pdb2-$chid2\n";
}
}

sub get_jail3{
my (@jail)=@_;
for(my $i=0;$i<$#jail;$i++){
	my ($uniprotAC1,$uniprotAC2,$shift1,$shift2);
	if($jail[$i] =~ /^#/){next;}# skip line
	if($jail[$i] =~ /^$/){next;}# skip line
	my @int = split(/\t/,$jail[$i]);
	my $pdb = $int[1];
	my $chid1 = $int[2];
	my $pdbch1 = $pdb.$chid1;
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
	printf OUT3 "$pdb\t$chid1\t$uniprotAC1\t$shift1\n";
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
