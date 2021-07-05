#!/usr/bin/perl -w
use strict;

my $est=shift;
my $freq_file = shift;
open(FREQ, $freq_file) or die "Cannot open $freq_file: $!";
my @content=<FREQ>;
close(FREQ);
chomp @content;

my %lib=();
foreach (@content){
	my ($chr, $transcript, $str)=(split/\s+/);
	my @snp=split/,/,$str;
	foreach (@snp){
		my ($maf,$pos,$fold)=split/_/;
		$pos++; # from 0 based to 1 based
		$lib{$transcript}->{$pos}=1;
	}
}

open(IN, $est) or die "Cannot open $est: $!";
my @rows=<IN>;
close(IN);
chomp @rows;
open(OUT, ">estSFS_full.txt");
foreach my $row (@rows){
	my ($transcript, $pos)=(split/\s+/,$row)[0,1];
	print OUT "$row\n" if(exists $lib{$transcript}->{$pos});
}
close(OUT);
