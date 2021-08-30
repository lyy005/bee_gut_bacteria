#!/usr/bin/perl -w
use strict;

die "perl $0 [cds_from_genomic.fna] [name list] [output file]\n" unless @ARGV == 3;

my %hash;

open (IN, "$ARGV[0]") or die "$ARGV[0] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
$/=">";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\n+/;
	my $id = shift @line;
	my $seq = join "", @line;
	
	if($id =~ /\[locus_tag=(\w+\_\w+)\]/){
		my $gene_id = $1;
		$hash{$gene_id} = $seq;
	}
}
$/ = "\n";
close IN;

open (LS, "$ARGV[1]") or die "$ARGV[0] $!\n";
<LS>;
while(<LS>){
	chomp;
	my @line = split;
	if($hash{$line[0]}){
		print OUT ">$line[0]\_$line[2]\n$hash{$line[0]}\n";
	}else{
		print "Not found: $line[0]\t$line[2]\n";
	}
}
close LS;
close OUT;
