#!/usr/bin/perl -w
use strict;

my $infile=shift;
my %map=();
open(IN, $infile) or die "Cannot open $infile: $!";
while(<IN>){
	chomp;
	my ($chr, $transcript, $freq_str)=(split/\s+/);
	my @freqs=split/,/,$freq_str;
	foreach (@freqs){
		my ($freq, $pos, $degenerate)=split/_/;
		push @{$map{$chr}->{$transcript}{$degenerate}},$freq if($degenerate eq '0' or $degenerate eq '4');
	}
}
close(IN);

open(OUT1, ">freq_0.txt");
open(OUT2, ">freq_4.txt");
print OUT1 "CHR\tTranscript\tFreq\n";
print OUT2 "CHR\tTranscript\tFreq\n";

foreach my $chr (sort keys %map){
	foreach my $transcript (sort keys %{$map{$chr}}){
		if(exists $map{$chr}->{$transcript}{'0'}){
			foreach (@{$map{$chr}->{$transcript}{'0'}}){
				print OUT1 "$chr\t$transcript\t",$_,"\n";
			}
		}
		if(exists $map{$chr}->{$transcript}{'4'}){
			foreach (@{$map{$chr}->{$transcript}{'4'}}){
				print OUT2 "$chr\t$transcript\t",$_,"\n";
			}
		}
	}
}