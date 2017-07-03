#!usr/bin/perl

#!/usr/bin/env python

# Copyright 2016, Rohan Dandage <rraadd_8@hotmail.com>
# This program is distributed under General Public License v. 3.  

use strict;
use warnings;
use Bio::DB::Sam;
use List::Util qw(sum);
use 5.014;
use POSIX qw/ceil/;
use Bio::Perl;
use Statistics::Basic qw(:all);
use File::Touch;
use File::Slurp;
use Bio::Seq;
use Bio::SeqIO;

 my $sam = Bio::DB::Sam->new(-fasta=>"$ARGV[0]", # ref 
                               -bam  =>"$ARGV[1]",	# sorted indexed bam
			       -expand_flags  => 1);

 my @alignments = $sam->get_features_by_location(-seq_id => 'gmr-wt');

my @snps; 
my $mutmatdatatmp1; #for mutmat
my @mutmatdatatmp2; #for mutmat
my @mutmatdata;

 for my $a (@alignments) {
    
    my $start  = $a->start;
    my $end    = $a->end;
    my $ref_dna= $a->dna;
    my $query_id   = $a->query->display_name;
    my $query_dna    = $a->query->dna;   
    my $cigar     = $a->cigar_str;
    my @scores    = $a->qscore;     # per-base quality scores
	my $tags = $a->aux;    
	my $basesincodon=3;
		@mutmatdatatmp2=""; #for mutmat	
	
# remove I D N
# if S > if S@3' 
# else normal
if ($tags !~ /NM:i:0/) {
if (($query_dna !~ /N/) && ($cigar !~ /D/) && ($cigar !~ /I/)){

	my @cig=split '',$cigar;
	my @ref_base=split '',$ref_dna;
	my @query_base=split '',$query_dna;
	my $S='S';	
	my $Sidx=index($cigar,$S);
	my @Snum;
#	print "$cigar\t$Sidx\n"
if (($cigar =~ /S/) && ($Sidx <= 3)) {
	my $s= substr($cigar,0,$Sidx);
	push (@Snum, $s);
}else {
	my $s=0;
	push (@Snum, $s);
}
			for (my $i=0;$i<@ref_base;$i++) {
			my $idx=$i+$start;
			if ($ref_base[$i] ne $query_base[$i+$Snum[0]]) {

					my $codonidxtemp=$idx/$basesincodon;
					my $codonidx = ceil($codonidxtemp);
					my $codonposi=$idx % $basesincodon; # 1 2 3(0)
								my @refcodon;
									if (($codonposi == 1) && (exists $ref_base[$i+1]) && (exists $ref_base[$i+2])) {				
									my $refcodontmp = join '',$ref_base[$i] , $ref_base[$i+1] , $ref_base[$i+2];
									push (@refcodon, $refcodontmp);
									
					}				elsif (($codonposi == 2) && (exists $ref_base[$i-1]) && (exists $ref_base[$i+1])) {
									my $refcodontmp = join '',$ref_base[$i-1] , $ref_base[$i] , $ref_base[$i+1];
									push (@refcodon, $refcodontmp);
					}				elsif (($codonposi == 0) && (exists $ref_base[$i-1]) && (exists $ref_base[$i-2])) {
									my $refcodontmp = join '',$ref_base[$i-2] , $ref_base[$i-1] , $ref_base[$i];
									push (@refcodon, $refcodontmp);
					}
								my @querycodon;
								my @scoresavg_overcodon=0; #scores are per query base for clipped bases also
									if (($codonposi == 1) && (exists $query_base[$i+1+$Snum[0]]) && (exists $query_base[$i+2+$Snum[0]])) {				
									my $querycodontmp = join '',$query_base[$i+$Snum[0]] , $query_base[$i+1+$Snum[0]] , $query_base[$i+2+$Snum[0]];
									push (@querycodon, $querycodontmp);
									my $scoresavg_overcodontmp=mean ($scores[$i+$Snum[0]] , $scores[$i+1+$Snum[0]] , $scores[$i+2+$Snum[0]]);
									@scoresavg_overcodon=();
									push (@scoresavg_overcodon,$scoresavg_overcodontmp);
									
					}				elsif (($codonposi == 2) && (exists $query_base[$i-1+$Snum[0]]) && (exists $query_base[$i+1+$Snum[0]])) {
									my $querycodontmp = join '',$query_base[$i-1+$Snum[0]] , $query_base[$i+$Snum[0]] , $query_base[$i+1+$Snum[0]];
									push (@querycodon, $querycodontmp);
									my $scoresavg_overcodontmp = mean( $scores[$i-1+$Snum[0]] , $scores[$i+$Snum[0]] , $scores[$i+1+$Snum[0]]);
									@scoresavg_overcodon=();
									push (@scoresavg_overcodon,$scoresavg_overcodontmp);

					}				elsif (($codonposi == 0) && (exists $query_base[$i-1+$Snum[0]]) && (exists $query_base[$i-2+$Snum[0]])) {
									my $querycodontmp = join '',$query_base[$i-2+$Snum[0]] , $query_base[$i-1+$Snum[0]] , $query_base[$i+$Snum[0]];
									push (@querycodon, $querycodontmp);
									my $scoresavg_overcodontmp = mean( $scores[$i-2+$Snum[0]] , $scores[$i-1+$Snum[0]] , $scores[$i+$Snum[0]]);
									@scoresavg_overcodon=();
									push (@scoresavg_overcodon,$scoresavg_overcodontmp);
					}
					my $refaatmp = translate(@refcodon);
					my $refaa      = $refaatmp->seq();
					my $queryaatmp = translate(@querycodon);
					my $queryaa      = $queryaatmp->seq();
					# synonymous mut
					my @syn;
					if ($refaa eq $queryaa) {
					my $s='y';
					push (@syn, $s);					
					} else {
					my $s='';
					push (@syn, $s);
					} 
					# aa to int
					my $aas='ACDEFGHIKLMNPQRSTVWY*'; 					
					my $queryaainttmp=index($aas,$queryaa);
					my $queryaaint=$queryaainttmp+1;
#
			if ($scoresavg_overcodon[0] >= 25) {
			$mutmatdatatmp1="$query_id\t$codonidx\t$queryaaint\t$refaa\t$queryaa\n";
			push (@mutmatdatatmp2, $mutmatdatatmp1);#, "\n";
}
}
} # for refbase
#print "@mutmatdatatmp2";
#---- mutmatdata ------ 
my $item;
my  %seen = ();
 foreach $item (@mutmatdatatmp2) {
   push(@mutmatdata, $item) unless $seen{$item}++; # dedup
}
@mutmatdatatmp2=();
#-----------------------
} # I D N
} # NM
} #alignment

#---- mutmat------
	my @mutmat;
	my $mutmatrixtmp;
#print "@mutmatdata";
foreach my $r (1 .. 21) {
foreach my $c (1 .. 178) {
my $mutpattern= "$c\t$r";
#print "$mutpattern";
my $count=()=grep(/\t$mutpattern\t/,@mutmatdata); 
$mutmatrixtmp="$count\t";
if ($c==178) {
push @mutmat, $mutmatrixtmp, "\n";
} else{ 
push(@mutmat, $mutmatrixtmp);
}
}
}
#---- gg stuff------
my $ggr=6;
my $ggcount=0;
my $ggitem_prev;
my $i;
my $ggitem;
my @ggmutmat="";

# wt G     12    20    21    25    33    48    52    68    89   115   121   138   149   161
for (my $ggc = 2; $ggc <= 178; $ggc += 2) {

# 1] wt G @ even
if (($ggc==12) || ($ggc==20) || ($ggc==48) || ($ggc==52) || ($ggc==68) || ($ggc==138)) {
			my 	$ggc_prev=$ggc-1;                        
			my	$pattern="$ggc_prev\t$ggr";
	$ggcount=()=grep(/\t$pattern\t/,@mutmatdata); 

# 2] wt G @ odd
} elsif 	(($ggc == 22) || ($ggc == 26) || ($ggc == 34) || ($ggc == 90) || ($ggc == 116) || ($ggc == 122)  || ($ggc == 150)  || ($ggc == 162)) {
			my	$pattern="$ggc\t$ggr";
	$ggcount=()=grep(/\t$pattern\t/,@mutmatdata); 
#print "$pattern\t$count\n";

# 3] "proper" GG's
} else {	for $i (0 .. $#mutmatdata) {
				$ggitem=$mutmatdata[$i];			
			my $ggc_prev;
			my $read_id;
			my $read_id_prev;
			
			my @ggitem_split;
			my @ggitem_prev_split;
			if ($ggitem =~ /\t$ggc\t$ggr\t/) {    # check in a line prev line
				@ggitem_split=split '\t',$ggitem;
				$read_id=$ggitem_split[0];         # get read_id 
				$ggc_prev=$ggc-1;                        
				if ($ggitem_prev =~ /\t$ggc_prev\t$ggr\t/) {     	#check pattern in prev line  
					@ggitem_prev_split=split '\t',$ggitem_prev;
					$read_id_prev=$ggitem_prev_split[0];          	# get read_id of prev line
						if ($read_id_prev eq $read_id) {  				# if the read id of both are same
							$ggcount++;                           			# count occasions.. 1 per occasion
						}
				}
			}
		$ggitem_prev=$ggitem;
		}
}
#print "$ggcount\n";
push (@ggmutmat, "$ggcount\n");
#empty 
$ggitem_prev='';
$ggcount=0;
}
#print "@ggmutmat";

# total 89 values
#------------------------------
#--------------normat----------

#------------------------------
#--------------export-----------
#my @files = ('mutmat','coverage');
#my $count = touch(@files);
  write_file( "$ARGV[1]mutmat", @mutmat ) ;
  write_file( "$ARGV[1]ggmutmatgg", @ggmutmat ) ;
#------------------------------
