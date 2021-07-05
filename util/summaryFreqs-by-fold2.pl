#!/usr/bin/perl -w
use strict;

my $freq_file=shift;
open(IN, $freq_file) or die "Cannot open $freq_file: $!";
my @content=<IN>;
chomp @content;
close(IN);

print STDOUT "Chr\tTranscript\tf0\tf4\tfN\tothers\n";
foreach (@content){
	my ($chr, $transcript, $freq_col)=split/\t/;
	print STDOUT "$chr\t$transcript\t";
	unless(defined $freq_col and $freq_col ne ''){
		print STDOUT "0\t0\t0\t0\n"; 
		next;
	}
	
	my @freqs=split/,/, $freq_col;
	my %fold_list=();
	my ($num_snp0, $num_snp4, $num_snpNull,$num_snp123)=(0,0,0,0);
	foreach (@freqs){
		my ($freq, $pos, $fold)=split/_/;
		# $freq=sprintf("%.2f", $freq);
		# if(exists $freq_list{$fold}->{$freq}){
		# 	$freq_list{$fold}->{$freq}++;
		# }else{
		# 	$freq_list{$fold}->{$freq}=1;
		# }
		if($fold == 0){
			$num_snp0++;
		}elsif($fold == 4){
			$num_snp4++;
		}elsif($fold == -1){
			$num_snpNull++;
		}else{
			$num_snp123++;
		}
	}
	print STDOUT "$num_snp0\t$num_snp4\t$num_snpNull\t$num_snp123\n";
}

# print  "\t",join("\t", sort {$a<=>$b} keys %{$freq_list{0}}),"\n0-fold:\t";
# foreach (sort {$a<=>$b} keys %{$freq_list{0}}){
# 	print $freq_list{0}->{$_},"\t";
# }
# print "\n4-fold:\t";
# foreach (sort {$a<=>$b} keys %{$freq_list{4}}){
# 	print $freq_list{4}->{$_},"\t";
# }
# print "\n";