#!/usr/bin/perl -w
use strict;

my $freq_file = shift;
my $anc_file = shift;
open(IN, $anc_file) or die "Cannot open $anc_file: $!";
$_=<IN>;
my @content =<IN>;
close(IN);
chomp @content;
my %map=();
foreach (@content){
	my ($transcript, $pos, $anc, $der_nb, $total_nb) = split/\t/;
	$pos--; # pos in frequency.txt is from 0 indexed
	my $der_freq = undef;
	if($anc ne 'NA'){
		if($total_nb == 0){
# 			print STDERR "$transcript, $pos, $anc, $der_nb, $total_nb\n";
		}else{
			$der_freq = $der_nb / $total_nb;
		}
	}
	
	$map{$transcript}->{$pos}= $der_freq;
}

open(FREQ, $freq_file) or die "Cannot open $freq_file: $!";
while(<FREQ>){
	chomp;
	my ($chr, $transcript, $info)=split/\t/;
	my @snp=split/,/,$info;
	my @str=();
	foreach (@snp){
		my ($freq, $pos, $fold)=split/_/;
		if(exists $map{$transcript}->{$pos} and defined $map{$transcript}->{$pos}){
			my $der_freq = $map{$transcript}->{$pos};
			if($der_freq > 0.5){
				$freq = $der_freq;
			}
		}else{
			$fold = -1;
		}
		push @str, $freq.'_'.$pos.'_'.$fold;
	}
	print "$chr\t$transcript\t",join("\t",@str),"\n";
}
close(FREQ);