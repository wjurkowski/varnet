#!/usr/bin/perl -w
use strict;
use warnings;

#get interfaces from 3did
my @threedid=open_file($ARGV[3]);
my $interactome="3did_interactome.tsv";
open(OUT, ">$interactome");
my $str="NO";
my $ni=0;
my (@int,@in,$nc,@cmap);
#parse interactions from 
#Other releases will include other sources of interactions
for(my $i=0;$i<$#threedid+1;$i++){
	my @int = split(/\t/,$threedid[$i]);
	my $id = $int[0];
	$id =~ s/\W//;	
	if($id eq "ID"){
		$ni++;
		$in[$ni][0]=$int[1];#name domain 1
		$in[$ni][1]=$int[2];#name domain 2
		$in[$ni][2]=substr($int[3],0,index($int[3],"."));#PFAM 1 id	
		$in[$ni][3]=substr($int[4],0,index($int[4],"."));#PFAM 2 id
		$str="NO";
	} 
	elsif($id eq "3D"){
		$in[$ni][4]=$int[1];#pdb id
		my @tab1=split(":",$int[2]);
		$in[$ni][5]=$tab1[0];#domain 1 chain id
		my @tab2=split("-",$tab1[1]);
		$in[$ni][6]=$tab2[0];#domain 1 start 
		$in[$ni][7]=$tab2[1];#domain 1 end
		my @tab3=split(":",$int[3]);
		$in[$ni][8]=$tab3[0];#domain 2 chain id
		my @tab4=split("-",$tab3[1]);
		$in[$ni][9]=$tab4[0];#domain 2 start
		$in[$ni][10]=$tab4[1];#domain 2 end
		$in[$ni][11]=$int[4];#score
		$in[$ni][12]=$int[5];#zscore
		$str="YES";		
	}	
	else{
	 if($str eq "YES"){
		$nc++;#count interacting AA pairs in given interface
		#pdb id; domain 1; PFAM id 1; chain id 1; domain s 1; domain e 1; domain 2; PFAM id 2; pdb id 2; chain id 2; domain s 2; domain e 2; AA 1; AA id 1; AA 2; AA id 2; score; zscore
		printf OUT "$in[$ni][4]\t$in[$ni][0]\t$in[$ni][2]\t$in[$ni][5]\t$in[$ni][6]\t$in[$ni][7]\t$in[$ni][1]\t$in[$ni][3]\t$in[$ni][8]\t$in[$ni][9]\t$in[$ni][10]\t$int[2]\t$int[0]\t$int[3]\t$int[4]\t$in[$ni][11]\t$in[$ni][12]\n";
	 }
	}
}