#!/usr/bin/perl -w
use strict;
use constant CUT => 0.8;
# use constant MIN_CUT => 0.2;

my $res_file = shift;
# my $chr_sel = shift;
my $coord_file = shift;
my $num_outgr = shift;
open(IN, $coord_file) or die "Cannot open $coord_file: $!";
my @coord = <IN>;
close(IN);
chomp @coord;
my @map=();
my $switch = 0;
foreach (@coord){
	my ($transcript, $pos, $ingr) = (split/\s+/)[0,1,2];
	my @allele_nb = split/,/,$ingr;
	my %table =(
			'A' => $allele_nb[0],
			'C' => $allele_nb[1],
			'G' => $allele_nb[2],
			'T' => $allele_nb[3]
		);
	
# 	my @major_alleles = sort{$table{$b} <=> $table{$a}} keys %table;
# 	my $major_allele = undef;
# 	if($table{$major_alleles[0]} == $table{$major_alleles[1]}){
# 		$major_allele = 'N'; ## equal counts
# 	}else{
# 		$major_allele = $major_alleles[0];
# 	}
	
	push @map,[$transcript, $pos, \%table]; 
}

open(IN, $res_file) or die "$res_file: $!";
my @content = <IN>;
close(IN);
chomp @content;
splice(@content,0,8);
die "Missing ancestral?" if($#content ne $#map);

my ($num_snp, $miss) = (0,0);
print STDOUT "transcript\tpos\tanc\tderived_allele_nb\ttotal_allele_nb\n";
foreach my $index (0..$#content){
	$num_snp++;
	my ($transcript, $pos, $ref_count) = @{$map[$index]};
	my @cols = split/\s+/,$content[$index];
	my $prob_major_anc = $cols[2];
	my $anc = undef;
	
	my $major_allele = (sort {$$ref_count{$b} <=> $$ref_count{$a}} keys %{$ref_count})[0]; # alleles with equal counts will fail in prob cutoff
	
	if($prob_major_anc >= CUT){
		$anc = $major_allele;
	}else{
		splice(@cols,0,3);
		foreach my $allele (('A','C','G','T')){
			my $prob = undef;
			if($num_outgr == 3){
				my @probs = splice(@cols,0,4);   ## for three outgrs
				$prob = $probs[0] + $probs[1] + $probs[2] + $probs[3]; 
			}elsif($num_outgr == 2){
				$prob = shift @cols # for two outgr
			}
			
			if($prob >= CUT){
				$anc = $allele;
				last;
			}   
		}
	}
	
	if(defined $anc){
		my $major_allele_nb = $$ref_count{$anc};
		my $total_allele_nb = $$ref_count{'A'} + $$ref_count{'C'} + $$ref_count{'G'} + $$ref_count{'T'};
		my $derived_allele_nb = $total_allele_nb - $major_allele_nb;
		print STDOUT "$transcript\t$pos\t$anc\t$derived_allele_nb\t$total_allele_nb\n";
	}else{
		my $total_allele_nb = $$ref_count{'A'} + $$ref_count{'C'} + $$ref_count{'G'} + $$ref_count{'T'};
		print STDOUT "$transcript\t$pos\tNA\tNA\t$total_allele_nb\n";
		$miss++;
	}	
}

print STDERR "$num_snp SNPs in total and $miss failed \n";








