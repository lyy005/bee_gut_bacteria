#!/usr/bin/perl -w
use strict;

die "perl $0 [S_wkB2_Amellifera.GC_geneID_powellID.list] [ortholog.table.col.sort] [output file]\n" unless @ARGV == 3;

my %hash;
my %spp;

open (IN, "$ARGV[1]") or die "$ARGV[1] $!\n";
open (OUT, ">$ARGV[2]") or die "$ARGV[2] $!\n";
<IN>;
while(<IN>){
	chomp;
	my $all = $_;
	my @line = split /\s+/;
	my @spp = split /\_/, $line[1];
	pop @spp;
	my $spp = join "_", @spp;
	
#print "$line[0]\t$spp\n";
	$hash{$line[0]}{$spp} ++;
	$spp{$spp} = 1;
}
close IN;

my %cate;
open (LS, "$ARGV[0]") or die "$ARGV[0] $!\n";
while(<LS>){
	chomp;
	my @line = split;

	if($cate{$line[-1]}){
		my @id = split /\_/, $line[0];
		pop @id;
		my $new_id = join "_", @id;
		$cate{$line[-1]} = $new_id."_".$cate{$line[-1]};
	}else{
		$cate{$line[-1]} = $line[0];
	}
#print "$line[-1]\t$line[0]\n";
}
close LS;

print OUT "GC\tfunction\tcomb\tspp\tnum\n";
foreach my $k (sort keys %cate){
	my $cate = $cate{$k};
	my @cate = split /\_/, $cate;
	my $cate_func = pop @cate;
	my $cate_ids = join "_", @cate;

	foreach my $l (sort keys %spp){
		if($hash{$k}{$l}){
			print OUT "$k\t$cate{$k}\t$cate_func\_$k\_$cate_ids\t$l\t$hash{$k}{$l}\n";
		}else{
			print OUT "$k\t$cate{$k}\t$cate_func\_$k\_$cate_ids\t$l\t0\n";
		}
	}
}
close OUT;
