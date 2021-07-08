#!/usr/bin/perl -w
use strict;
use File::Basename;
use File::Find;
use Cwd qw{abs_path};
use constant GT_COL=> 9; # start position of individual info
use constant NLINE => 100000;
my $gff=shift;
my $indir=shift;
my $outdir=shift;

my %pos_map=();

my $gff_fh;
if($gff=~/\.gz$/){
	open($gff_fh, "zcat $gff |") or die "Cannot open $gff: $!";
}else{
	open($gff_fh, $gff) or die "Cannot open $gff: $!";
}

print STDERR "Reading GFF\n";
while(<$gff_fh>){
	chomp;
	next if(/^#/);
	my ($chr, $feat, $start, $end)= (split/\s+/)[0,2,3,4];
	next unless($feat eq 'CDS');
	foreach ($start .. $end){
		$pos_map{$chr}->{$_}=1;
	}
}
close($gff_fh);
#
print STDERR "Start to convet VCF\n";
my @filenames=get_filenames_recursive($indir, 'vcf.gz');
foreach my $file (@filenames){
	print STDERR "\nProcessing $file ... \n";
	open(my $fh, "gzip -d -c $file|") or die "Cannot open  $file: $!";
	## header
	my @ind_ids=();
	while(<$fh>){
		next if(/^#{2,}/);
		if(/^#CHROM/){
			chomp;
			my $header_line=$_;
			my @header_cols=split/\t/,$header_line;
			@ind_ids=splice(@header_cols, GT_COL);
			last;
		}
	}
	my @info_lines = map~~<$fh>, 1..NLINE;
	
	while($info_lines[0]){
		parse_vcf(\@info_lines, \@ind_ids);
		@info_lines = map~~<$fh>, 1..NLINE;
		# print STDERR "Done!";
	}
	close($fh);
}
print STDERR "Done!\n";
	
sub parse_vcf{
	my ($ref_lines, $ind_ids_ref)=@_;
	my @snp_list=();
	my @ind_ids=@$ind_ids_ref;
	# print STDERR "parsing";
	OUT1:foreach (@$ref_lines){
		last if($_ eq ''); # after EOF empty assign to array
		chomp;
		my @cols=split/\t/;
		my ($chr_id, $pos, $ref, $alt, $info, $format)=@cols[0,1,3,4,7,8];
		next unless(exists $pos_map{$chr_id}->{$pos}); # keep only the CDS regions
		
		next if($alt =~/\*/ || $info=~/^INDEL/ || length($ref) > 1);# alt allele is missing due to a upstream deletion || INDEL in either ref or alt. 
		my @alleles=();
		if($alt eq '.'){
			$alt = $ref;
			push @alleles, $ref;
		}else{
			my @alt_alleles=split/,/,$alt;
			foreach (@alt_alleles){
				next OUT1 if(length($_)>1); 
			}
			push @alleles, ($ref, @alt_alleles);
		}
		
		
		my @inds=splice(@cols,GT_COL);				
		die "Missing columns in \#snp$chr_id, $pos" unless($#ind_ids == $#inds);
		# print "here\n";
		foreach my $ind_index (0..$#ind_ids){
			my @genotype_allele_indice=split/\/|\|/,(split/:/,$inds[$ind_index])[0]; # remove the AD, and split GT by / (unphased) or | (phased)
			my $genotype=determine_genotype(\@alleles, \@genotype_allele_indice);
			if(defined $genotype){
				$genotype='N' if($genotype eq '*');
			}else{
				$genotype='N';
			}
			next if($genotype eq $ref);
			my $print_snp_str=$chr_id."\t".$pos."\t".$ref."\t".$genotype;
			push  @{$snp_list[$ind_index]}, $print_snp_str;
		}
	}
	print STDERR ".";
	foreach my $ind_index (0..$#snp_list){
		next unless(defined $snp_list[$ind_index]);
		my $output=$outdir.'/'.$ind_ids[$ind_index].'.snp';
		open(OUT, '>>', $output) or die "Cannot write into $output: $!";
		print OUT join("\n", @{$snp_list[$ind_index]}),"\n";				
		close(OUT);
	}
}

sub extract_field_value{
	my $str=shift;
	my $field=shift;
	my ($value)=$str=~/$field=(-*\d*\.*\d+e*[\-|+]*\d*)/i;
	if(defined $value){
		return($value);
	}else{
		return('NA');
	}
}

sub determine_genotype{
	my ($alleles_arr, $index_arr)=@_;
	my $gt='';
	foreach (@$index_arr){
		my $allele='';
		if(/\d+/){
			$allele=$$alleles_arr[$_];
		}elsif($_ eq '.'){ # missing, then be 'N'
			$allele='N';
		}else{
			return(undef); ## cannot be anything other than number or '.'
		}

		if($allele=~ /INS/){  # insert relative to ref
			$allele='?'; # '+' not allowed in bioperl so we use ? for insertion
		}elsif($allele=~/DEL/){ # deletion relative to ref
			$allele='-';
		}

		if($gt eq $allele){
			# do nothing
		}else{
			$gt.=$allele;
		}
	}
	# $gt=uc($gt);
	my %IUPAC_DNA=('AG' => 'R', 'GA' => 'R',
				   'CT' => 'Y', 'TC' => 'Y',
				   'GC' => 'S', 'CG' => 'S',
				   'AT' => 'W', 'TA' => 'W',
				   'GT' => 'K', 'TG' => 'K',
				   'AC' => 'M', 'CA' => 'M',
				   'CGT' => 'B','CTG' => 'B','GTC' => 'B', 'GCT' => 'B', 'TCG' => 'B', 'TGC' => 'B',
				   'AGT' => 'D','ATG' => 'D','GTA' => 'D', 'GAT' => 'D', 'TAG' => 'D', 'TGA' => 'D',
				   'ACT' => 'H','ATC' => 'H','CTA' => 'H', 'CAT' => 'H', 'TAC' => 'H', 'TCA' => 'H',
				   'ACG' => 'V','AGC' => 'V','CGA' => 'V', 'CAG' => 'V', 'GAC' => 'V', 'GCA' => 'V',
				    );
	if($gt eq ''){
		$gt=undef;
	}elsif(length($gt)==1){
		# do nothing
	}elsif(exists $IUPAC_DNA{$gt}){
		$gt=$IUPAC_DNA{$gt};
	}elsif($gt=~/[ATGC]{4}/){
		$gt='N';
	}else{
		$gt=undef;
	}
	return($gt);
}

sub get_field_index_by_format{
	my $format=shift;
	my @fields=split/:/,$format;
	my ($ad_index, $dp_index, $geno_qual_index, $pl_index)=(undef, undef, undef, undef);
	my $index=0;
	foreach (@fields){
		$ad_index=$index if($_ eq 'AD');
		$dp_index=$index if($_ eq 'DP');
		$geno_qual_index=$index if($_ eq 'GQ' || $_ eq 'RGQ');
		# $pl_index=$index if($_ eq 'PL');
		$index++;
	}
	# if(defined $ad_index and defined $dp_index and defined $geno_qual_index){
		return(($ad_index, $dp_index, $geno_qual_index));
	# }else{
		# die "Unknown FORMAT fileds: $format";
	# }
}

sub get_filenames_recursive{
	my ($parent_path, $filter)=@_;
	$parent_path=abs_path($parent_path);
	my @filenames=();
	find sub{
		push @filenames, $File::Find::name if(-f $File::Find::name && grep(/\.$filter$/i, $_));
	}, $parent_path;
	return(@filenames);
}
