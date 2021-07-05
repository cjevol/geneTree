#!/usr/bin/perl -w
use strict;

my $infile =shift;
open(IN, $infile) or die "Cannot open $infile: $!";
my @content=<IN>;
close(IN);
chomp @content;

print STDOUT "transcript\tpos\tanc\tderived_allele_nb\ttotal_allele_nb\n";
foreach (@content){
	my ($transcript, $pos, $ingr, $outgr) = (split/\t/)[0,1,2,3];
	
	my @allele_nb = split/,/,$ingr;
	my %ingr_counts =(
		'A' => $allele_nb[0],
		'C' => $allele_nb[1],
		'G' => $allele_nb[2],
		'T' => $allele_nb[3]
	);
	
	@allele_nb = split/,/,$outgr;
	my %outgr_counts =(
		'A' => $allele_nb[0],
		'C' => $allele_nb[1],
		'G' => $allele_nb[2],
		'T' => $allele_nb[3]
	);

# 	my @in_alleles_srt = sort {$ingr_counts{$b} <=> $ingr_counts{$a}} keys %ingr_counts;	
	my @out_alleles_srt = sort {$outgr_counts{$b} <=> $outgr_counts{$a}} keys %outgr_counts;	
	
	my $outgr_allele_fix = $out_alleles_srt[0];
	my $total_allele_nb = $ingr_counts{'A'} + $ingr_counts{'C'} + $ingr_counts{'G'} + $ingr_counts{'T'};
	my $anc_allele_nb = $ingr_counts{$outgr_allele_fix};
	my $derived_allele_nb = $total_allele_nb - $anc_allele_nb;
	print STDOUT "$transcript\t$pos\t$outgr_allele_fix\t$derived_allele_nb\t$total_allele_nb\n";
}

