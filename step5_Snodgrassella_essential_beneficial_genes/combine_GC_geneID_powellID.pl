#!/usr/bin/perl -w
use strict;

die "perl $0 [GC2geneID] [blast filter out geneID2powellID] [output file]\n" unless @ARGV == 3;

my %hash;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\s+/;
	
	$hash{$line[1]} = $line[0];
}
close IN;

open (LS, "$ARGV[1]") or die "$ARGV[0] $!\n";
while(<LS>){
	chomp;
	my @line = split;
	
	print OUT "$line[3]\t$line[2]\t$hash{$line[2]}\n";
}
close LS;
close OUT;
