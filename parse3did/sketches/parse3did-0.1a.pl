#!/usr/bin/perl -w
use strict;
use warnings;

#Files needed on input: 3did interactome from 3did_flat
if ($#ARGV != 0) {die "Program used with parameters [interactome]\n";}



#get interfaces from 3did
my @threedid=open_file($ARGV[0]);
my $interactome="3did_interactome.tsv";
open(OUT, ">$interactome");
my $str="NO";
my $ni=0;
my (@in);
#parse interactions from 
#Other releases will include other sources of interactions
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
	} 
	elsif($id eq "3D"){
		$ni++;
		$in[4]=$int[1];#pdb id
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
#print "hola lola $id\n";
	}	
	else{
	 if($str eq "YES"){
		#print interacting AA pairs in given interface. 
		#Interface is defined as instance of pfam domain pairs in contact in particular PDB file. Given PDB file may include multiple interfaces representing same or different pfam-pfam interaction
		#Interface # is to identify same pfam-pfam in unique pdb
		#pdb id; interface #; domain 1; PFAM id 1; chain id 1; domain s 1; domain e 1; domain 2; PFAM id 2; pdb id 2; chain id 2; domain s 2; domain e 2; AA 1; AA id 1; AA 2; AA id 2; score; zscore
		printf OUT "$in[4]\t$ni\t$in[0]\t$in[2]\t$in[5]\t$in[6]\t$in[7]\t$in[1]\t$in[3]\t$in[8]\t$in[9]\t$in[10]\t$int[2]\t$int[0]\t$int[3]\t$int[4]\t$in[11]\t$in[12]\n";
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
