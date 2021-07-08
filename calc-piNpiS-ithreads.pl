#!/opt/local/bin/perl -w
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use lib './'; # 
use GenomeAnalysis;
use File::Basename;
use File::Path qw{remove_tree};
use Time::HiRes;
use Getopt::Long qw(GetOptions);
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(sum);
use Cwd;
use Benchmark;
use Pod::Usage qw(pod2usage);
use constant GT_COL => 9; ## NOTE CHANGE TO 8 FOR STANDARD VCF....
# use constant NUMH => 20; # number of header of VCF


################## initiate variables
my $num_seqs=0; # number of genes to pick randomly
my $output_align=0; # whether to print seq alignments
my $sample_size=10; # number of individual genomes to sample
my $tmp_workdir=getcwd().'/tmp'; # tmp working folder
my $num_threads=4;
my $snp_type='snp'; ### snp or vcf 
my $chr_sel=0;
my $ploid=1;

my $help=0;
if ((@ARGV==0) && (-t STDIN)) {
	 pod2usage(-verbose => 0, -message =>  "$0: Argument required.\n");
}
my $command_line=$0.' '.join(" ",@ARGV);

&GetOptions('chr=s'		=>	\$chr_sel, # chr# to cal, 0 for all
			'dip=i'     =>  \$ploid, # 1 or 2 accepted
            'format=s'   => \$snp_type,
            'gff=s'		=>	\my $gff_file, # gff
            'num=s' 	=> \$num_seqs, # random select and output num alignments to output (for PALM) 
            'out=s' 	=> \$output_align,
            'poly=s'    =>  \my $snpdir, # path to the snp
            'ref=s'		=>	\my $ref_genome, # ref genome
            'size=i'    => \$sample_size, # num genomes
            'threads=i'   => \$num_threads, 
            'work=s'		=>	\$tmp_workdir, # temp working directory
            'help|?'	=>	\$help) or pod2usage(-msg => 'Wrong options', -verbose => 1);

pod2usage(-verbose=>2) if $help;		
pod2usage(-msg => 'Please give a name for gff file', -verbose=>2) unless($gff_file);
pod2usage(-msg => 'Please give a name for genome fasta file', -verbose=>2) unless($ref_genome);
pod2usage(-msg => 'Please give a name for snp  directory', -verbose=>2) unless($snpdir);


print STDERR "\n####################################################################################################\n";
printf STDERR "%40s%20s\n","ANALYSIS","STARTED!";
print STDERR "####################################################################################################\n\n";
print STDERR "$command_line\n\n";

my $start_run=new Benchmark;
##########################################################################################################################################################
# 0 prepare the working directory
##########################################################################################################################################################
mkdir $tmp_workdir or die "Cannot build working directory $tmp_workdir: $!" unless(-d $tmp_workdir);

my ($gff_name)=fileparse($gff_file, qw{.GFF .gff .gff3 .gtf .gff.gz .gff3.gz .gtf.gz});
my ($ref_name)=fileparse($ref_genome, qw{.fasta .fa .fas .fasta.gz .fa.gz .fas.gz});
# $tmp_workdir.='/'.$ref_name.'-'.$gff_name;
mkdir $tmp_workdir or die "Cannot build working directory $tmp_workdir: $!" unless(-d $tmp_workdir);

my $snp_tmp_dir= $tmp_workdir.'/SNP';
mkdir $snp_tmp_dir or die "Cannot build folder $snp_tmp_dir: $!" unless(-d $snp_tmp_dir);

my $output_align_dir=undef;
if($output_align){
	$output_align_dir=$tmp_workdir.'/align';
	if(-d $output_align_dir){
		remove_tree($output_align_dir, {keep_root => 1}) || print STDERR "No files removed in $output_align_dir\n";
	}else{
		mkdir $output_align_dir || die "Cannot build $output_align_dir: $!";
	}
}
##########################################################################################################################################################
# step 1: read gff annotation file
##########################################################################################################################################################

my (%cds_feature_coord, %useful_coord_ref)=(); # {chr}->{trans_id}{strand}=[start, end]
$chr_sel =uc($chr_sel);
GenomeAnalysis::gff_parser($gff_file, $chr_sel, \%cds_feature_coord);

##########################################################################################################################################################
# step 1.2: sample genes
##########################################################################################################################################################

if($num_seqs){
	my ($counter,$dump_counter)=0;
	foreach my $chr (keys %cds_feature_coord){
		my @transcript_id=keys %{$cds_feature_coord{$chr}};
		$counter+=scalar @transcript_id;
		my @sampled_index=GenomeAnalysis::sample_array($num_seqs, 0..$#transcript_id);
		delete @transcript_id[@sampled_index]; # now it is only the dumped ids stored
		@transcript_id=grep($_, @transcript_id); # remove undef
		$dump_counter+=scalar @transcript_id;
		delete @{$cds_feature_coord{$chr}}{@transcript_id}; # delete the dumped 
		# print join("\n", keys %{$cds_feature_coord{$chr}}),"\n";
	}
	if($chr_sel){
		print STDERR "$counter transcripts have been identified in Chr$chr_sel in GFF\n";
		print STDERR $counter - $dump_counter, " transcripts have been sampled\n";
	}else{
		print STDERR "$counter transcripts have been identified in GFF\n";
		print STDERR $counter - $dump_counter, " transcripts have been sampled ($num_seqs for each chrom)\n";
	}
}else{
	my $counter=0;
	foreach my $chr (keys %cds_feature_coord){
		$counter+=scalar keys %{$cds_feature_coord{$chr}};
	}
	if($chr_sel){
		print STDERR "$counter transcripts have been identified in Chr$chr_sel in GFF\n";
	}else{
		print STDERR "$counter transcripts have been identified in GFF\n";
	}	
}

##########################################################################################################################################################
# step 2: sample and split snp files by chromosome
##########################################################################################################################################################

my @chr_ids=keys %cds_feature_coord;

my @snp_filenames=GenomeAnalysis::get_filenames_recursive($snp_tmp_dir, $snp_type); # check if snp files have been manually selected
if(@snp_filenames){
	print STDERR scalar @snp_filenames, " snp files have been found\n";
}else{
	### randomly pick individual genomes and split snp files by chromosomes
	if($snp_type eq 'snp'){
		@snp_filenames=GenomeAnalysis::split_snp_file($snpdir, $snp_type, $sample_size,$snp_tmp_dir, \@chr_ids, $num_threads);
	}elsif($snp_type =~/vcf/i){
		@snp_filenames=GenomeAnalysis::split_vcf_file($snpdir, $sample_size,$snp_tmp_dir, \@chr_ids, $num_threads);
	}else{

	}
	print STDERR "There are ",scalar @snp_filenames, " snp files after split\n";
}

##########################################################################################################################################################
# step 3. for each chr, get CDS, get SNP, then Pi
##########################################################################################################################################################
print STDOUT "Chr\tTranscript\tLength\tL0\tL4\tSNP0\tSNP4\tPi0\tPi4\n";
my $ref_fh;
if($ref_genome=~/\.gz$/){
	open($ref_fh, "zcat $ref_genome |") or die "Cannot open $ref_genome: $!";
}else{
	open($ref_fh, $ref_genome) or die "Cannot open $ref_genome: $!";
}
my $in=Bio::SeqIO->new(-fh => $ref_fh, -format => 'fasta');
# my @hash_arr=(); # contain shared hashes [shash{ind1=>pos=alt}, shash{ind2=>pos=alt}, .... ]
my @snp_list;
my @chr_results2print:shared;
my @n_ratio2print:shared;
my @freq_str2print:shared;

open(MYERR, ">missing_ratio.txt") or die "Cannot open fh for writing: $!";
print MYERR "Chr\tTranscript\tNum_N\tNum_SNP\n";
open(FREQ, '>frequency.txt') or die "Cannot open fh for writing: $!";
while(my $seqObj=$in->next_seq){
	my $chr_id=uc($seqObj->display_id);
	
	# if($chr_id=~/(\d+)/ or $chr_id=~/([XY])/i){ #only chromosomes with numbers and x, y
	# 	$chr_id=$1;
		next unless(!$chr_sel || $chr_id eq $chr_sel);
	# }else{
# 		next;
# 	}
	
	next unless(exists $cds_feature_coord{$chr_id});
	print STDERR "Processing Chr$chr_id with ", scalar keys %{$cds_feature_coord{$chr_id}}," transcripts\n";
	
	####################################################################################
	# 	# step 3.1 read snps
	#####################################################################################
	@snp_list=();
	@chr_results2print=();
	@n_ratio2print=();
	@freq_str2print=();
	
	my @snp_filenames_by_chr=grep(/\/$chr_id\-[^\/]+\.$snp_type$/i, @snp_filenames);
	next unless(@snp_filenames_by_chr);
	print STDERR "\tSNP files found: \n\t", join("\n\t", @snp_filenames_by_chr),"\n\n";
	if($snp_type eq 'snp'){
		foreach my $ind_index (0..$#snp_filenames_by_chr){
			my $snp_file_chr=$snp_filenames_by_chr[$ind_index];
			my ($ind_id)=fileparse($snp_file_chr, qw{.snp});
			# $snp_list[$ind_index]=shared_clone({});
			
			print STDERR "\tReading $snp_file_chr...";
			open(my $fh, '<', $snp_file_chr) || die "Cannot open $snp_file_chr: $!";
			my @content=<$fh>;
			chomp @content;
			close($fh);
			# GenomeAnalysis::process_array_ithread(\@content, $num_threads, \&process_snp_file, $ind_index);
			$snp_list[$ind_index]={};
			process_snp_file(\@content, $snp_list[$ind_index]);
			print STDERR "\tDone!\n\t", scalar keys %{$snp_list[$ind_index]}," snps were found in $ind_id\n";
		}
	}elsif($snp_type=~/vcf/i){
		die "One VCF file with multiple individuals is allowed!" if ($#snp_filenames_by_chr>0);
		my $snp_file_chr=$snp_filenames_by_chr[0];
		
		print STDERR "\tReading $snp_file_chr...";
		open(my $fh, '<', $snp_file_chr) or die "Cannot write into $snp_file_chr: $!";
		## header
		# my @info_lines = map~~<$fh>, 1..NUMH;
		my $header_line=<$fh>;
		chomp $header_line;
		my @header_cols=split/\t/,$header_line;
		my @ind_id=splice(@header_cols, GT_COL);
		push @snp_list, shared_clone({}) for @ind_id;
		
		my @content=<$fh>;
		chomp @content;
		close($fh);
		GenomeAnalysis::process_array_ithread(\@content, $num_threads, \&process_vcf_file);
		
		print STDERR "\tDone!\n\t", scalar @snp_list, " inds with ",scalar keys %{$snp_list[0]}, " snps have been found in $snp_file_chr\n";
		
	}else{
		#### need to add other format for snp files, e.g. tubular 
	}
	
	####################################################################################
	# 	# step 3.2 extract CDS and calculate PI
	#####################################################################################
	print STDERR "\tSTART to calculate nucleotide diversity...\n";
	my @transcript_id = keys %{$cds_feature_coord{$chr_id}};
	GenomeAnalysis::process_array_ithread(\@transcript_id, $num_threads, \&calculate_pi_per_cds, ($seqObj, $cds_feature_coord{$chr_id}, \@snp_list, $ploid, $output_align_dir,\@snp_filenames_by_chr));
	print STDERR "\tDone!\n";
	
	####################################################################################
	# 	# step 3.3 print results to file
	#####################################################################################
	print STDERR "\tPrint results to output!\n";
	for (@chr_results2print){
		print STDOUT "$chr_id\t", $_,"\n";
	}
	for (@n_ratio2print){
		print MYERR "$chr_id\t", $_,"\n";
	}
	for (@freq_str2print){
		print FREQ "$chr_id\t", $_,"\n";
	}	
	last if($chr_id eq $chr_sel);
}
close(MYERR);
close(FREQ);
my $end_run=new Benchmark;
print STDERR "\nJob took ", timestr(timediff($end_run,$start_run)),"\n";
print STDERR "\n####################################################################################################\n";
printf STDERR "%40s%20s\n","ANALYSIS","COMPLETED!";
print STDERR "####################################################################################################\n\n";


####################################################################################################################################################################################

sub calculate_pi_per_cds{
	my ($q, $seqObj, $cds_coords, $snp_list, $ploid, $output_align_dir, $snp_file_name_ref)=@_;
	my $out;
		
	while(defined(my $transcript_id=$q->dequeue())){
		my (%bad,%snp_sites, %snp_classified)=(); # to keep bad snp (missing 'N', or ambiguous 'R', 'W' in haploid data)
		my ($cds_seqObj, $cds_genome_coords_arr)=GenomeAnalysis::extract_CDS_by_coord($seqObj, $$cds_coords{$transcript_id}, $transcript_id);
		next unless(defined $cds_seqObj);
		my $cds_sequence=uc($cds_seqObj->seq);
		my $cds_length=$cds_seqObj->length;
		
		if(defined $output_align_dir){
			$out=Bio::SeqIO->new(-file => '>'.$output_align_dir.'/'.$transcript_id.'.fas', -format=>'fasta');
		}
		# update individual info
		
		my (%polymorphism, %degeneracy_fold, %ref_list, @freq)=();
		foreach my $ind (0..$#$snp_list){
			my $seq_name=basename($$snp_file_name_ref[$ind],'.snp');
			$seq_name=(split/-/,$seq_name)[1];

			my $ind_seq;
			{	# for safe lock the shared snp_list 
				# lock($$snp_list[$ind]);
				$ind_seq=GenomeAnalysis::get_indSeq_by_refSeq_polym($cds_sequence, $cds_genome_coords_arr, $$snp_list[$ind], \%polymorphism, \%ref_list, $ploid);
			}
			my $ind_seqObj=Bio::Seq->new(-id=>$seq_name, -seq=>$ind_seq, -alphabet=>'dna');
			
			if($output_align){
				$out->write_seq($ind_seqObj);
			}
			my @degeneracy_fold=GenomeAnalysis::get_degenerate_fold_by_seq($ind_seqObj);
			unless($cds_length == ($#degeneracy_fold +1) ){
				print STDERR "we have something wrong when getting degeneracy for individuals ", $cds_seqObj->id, " $ind\n";
				die "$!";
			}
			foreach (0..$#degeneracy_fold){
				$degeneracy_fold{$_}->{$degeneracy_fold[$_]}=1;
			}
		}
		# calculate pi and degeneracy for each position
		my (%degenerate_length, %degenerate_pi)=();
		my $actual_sample_size=($#$snp_list+1)*$ploid;
		
		# print STDERR "$transcript_id\t"; ####### checking
		foreach my $pos (0..($cds_length-1)){
			my @base_fold=keys %{$degeneracy_fold{$pos}};
			# $snp_sites{$pos}=1 if(exists $polymorphism{$pos});
			
			if(exists $ref_list{$pos}){ # make sure we only have alternative alleles here
				my $ref_allele=$ref_list{$pos};
				# my @alleles=keys %{$polymorphism{$pos}};####### checking
				# print STDERR $pos," ",$ref_allele, ":",join(",", @alleles),'['; ####### checking
				if(exists $polymorphism{$pos}->{$ref_allele}){
					if($ploid==1){# for haploid data, there should not be ref allele in alt alleles, i.e. no ambiguous alleles like R, W, K, M... If there is, we don't count this site
						$bad{$pos}=1;
						$snp_sites{$pos}=1;
						next;
					}
					delete $polymorphism{$pos}->{$ref_allele};
				}
			}
			
			### determine the fold
			die "We have trouble to get fold info in $transcript_id pos$pos" unless(@base_fold);
			my ($degerate_base, $pi_base)=();
			if($#base_fold>0){
				$degerate_base=-1;
			}else{
				$degerate_base=$base_fold[0];
			}
			
			# count Ln, Ls
			
			my @alleles=keys %{$polymorphism{$pos}}; # alt alleles
			if($#alleles>0 && $degerate_base!=-1){
				$bad{$pos}=1;
				$snp_sites{$pos}=1;
				next;
			} # if we have more than 1 alternative alleles in the position, ignore the sites even if we can determine the fold 
			
			if(exists $degenerate_length{$degerate_base}){
				$degenerate_length{$degerate_base}++;
			}else{
				$degenerate_length{$degerate_base}=1;
			}
			
			next unless($#alleles==0); # no polymorphism
			
			## calculate pi
			my $num_alt_allele=$polymorphism{$pos}->{$alleles[0]};
			my $num_ref_allele=$actual_sample_size-$num_alt_allele;
			my $base_num_diff=$num_alt_allele*$num_ref_allele;
			if($base_num_diff ==0){
				next; # all alleles are alt alleles ... so polymorphic sites sampled even though in all samples might exist
			}
			$snp_sites{$pos}=1;
			my $freq=$num_alt_allele/$actual_sample_size;
			if($freq>0.5){
				push @freq, (1-$freq).'_'.$pos.'_'.$degerate_base;
			}else{
				push @freq, $freq.'_'.$pos.'_'.$degerate_base;
			}
			# print STDERR $num_alt_allele,":",$num_ref_allele,"];"; ####### checking
			if(exists $degenerate_pi{$degerate_base}){
				$degenerate_pi{$degerate_base}+=$base_num_diff;
				$snp_classified{$degerate_base}+=1;
			}else{
				$degenerate_pi{$degerate_base}=$base_num_diff;
				$snp_classified{$degerate_base}=1;
			}
		}
		# print STDERR "\n"; ####### checking
		my $total_pairs=$actual_sample_size*($actual_sample_size-1)/2;
		my $res_str=$transcript_id."\t".$cds_length."\t";
		my $freq_str=$transcript_id."\t".join(",",@freq);
		# add Ln
		if(exists $degenerate_length{'0'}){
			$res_str.=$degenerate_length{'0'}."\t";
		}else{
			$res_str.="0\t";
		}
		# add Ls; for Nei 1986, S=S4+1/3*S2+2/3*S3; and N=n-S
		if(exists $degenerate_length{'4'}){
			$res_str.=$degenerate_length{'4'}."\t";
		}else{
			$res_str.="0\t";
		}
		# add snp
		if(exists $snp_classified{'0'}){
			$res_str.=$snp_classified{'0'}."\t";
		}else{
			$res_str.="0\t";
		}
		if(exists $snp_classified{'4'}){
			$res_str.=$snp_classified{'4'}."\t";
		}else{
			$res_str.="0\t";
		}
		## add pi_N
		if(exists $degenerate_pi{'0'}){
			$res_str.=$degenerate_pi{'0'}/$total_pairs."\t";
		}else{
			$res_str.="0\t";
		}
		## add pi_S
		if(exists $degenerate_pi{'4'}){
			$res_str.=$degenerate_pi{'4'}/$total_pairs."\t";
		}else{
			$res_str.="0\t";
		}
		## 
		my $num_bad_sites=scalar keys %bad;
		my $num_snp_sites=scalar keys %snp_sites;
		my $n_ratio_str="$transcript_id\t$num_bad_sites\t$num_snp_sites";
		
		{
			lock(@chr_results2print);
			push @chr_results2print, $res_str;
			push @n_ratio2print, $n_ratio_str;
			push @freq_str2print, $freq_str;
		}
	}
}

sub process_snp_file{
	my ($arr, $hash)=@_;
	foreach my $line (@$arr){
		my ($pos,$ref, $alt)=(split/\t/,$line)[1,2, 3];
		$$hash{$pos}=[$ref,$alt];
	}
}

sub process_vcf_file{
	my ($q)=@_;
	while(defined(my $line = $q->dequeue())){
		my @cols=split/\t/,$line;
		my ($pos, $ref, $alt)=@cols[1,3,4];
		next if($alt eq '*');# alt allele is missing due to a upstream deletion. 
		my @inds=splice(@cols,GT_COL);
		# next if($alt=~/DEL|INS|DUP|INV|CNV/); # we can remove this later
		my @alt_alleles=split/,/,$alt;
		my @alleles=();
		push @alleles, ($ref, @alt_alleles);
		foreach my $ind_index (0..$#inds){
			my @genotype_allele_indice=split/\/|\|/,(split/:/, $inds[$ind_index])[0]; # remove the AD, and split GT by / (unphased) or | (phased)
			my $genotype=GenomeAnalysis::determine_genotype(\@alleles, \@genotype_allele_indice);
			next unless(defined $genotype);
			lock($snp_list[$ind_index]);
			$snp_list[$ind_index]{$pos}=shared_clone([$ref, $genotype]);
		}
	}
}



__END__

###################################### usage
=head1 NAME

calc-piNpiS-ithreads.pl - This script calculates nucleotide diversity in 0-folded (Pi_N) and 4-folded sites (Pi_S) 

=head1 SYNOPSIS

 perl calc-piNpiS-ithreads.pl
           -chr <chr# for analysis only>
           -ref <ref genome> 
           -gff <GFF file> 
           -dip <ploid, default:1>
           -poly <folder contain polymorphism data [snp], can be recursive>
           -tmp <tmp working dir; default tmp> 
           -size <number of genomes to select by random; default 10>
           -num <number of genes to sample 0:all>
           -out <1: print sequence alignments; 0: Do not>
           -threads <number of threads>
           -format <snp | vcf |> we now accept .snp format [#chr\tPOS\tREF\tALT ...] or standard VCF
           -help show detail help message

=head1 OPTIONS

=over 4

=item B<--chr>

Chr# that chosen for analysis. Can be any number or X,Y; default:0 [all] 

=item B<--ref>

Reference genomic DNA sequences in 'FASTA' format. Mandatory parameter

=item B<--gff>

An annotation file in 'GFF3' format specifying 'CDS' regions

=item B<--poly>

The fold contains polymorphism [SNP] data, can contain recursive directory tree

=item B<--format>

The suffix of polymorphism files specifying the format [CASESENSITIVE]:
'.snp'-- Chr#[tab]Position[tab]Ref_allele[tab]Alt_allele[tab]Any number of fields for qualities  e.g. 1001 AT genome project  
'.vcf'-- Standard VCF file. One and Only one file contains SNP info for all individuals

=item B<--dip>

The number of ploid. Only 1 [default] or 2 accepted now
  
=item B<--size>

The size of subset of individual genomes to sample by random without replacement. So the number should not be larger than total number of individuals available
default: 10

=item B<--num>

The size of subset of CDS sequences to sample by random without replacement. 
default: 0 [all]

=item B<--out>

Whether to print out the CDS sequence alignments of sample individuals
if 1 [YES] print to folder 'align' under working directory;  default: 0 [NO] 


=item B<--tmp>

The working directory to store all files from intermediate steps which can be reused 
default: 'tmp' under current directory

=item B<--threads>

The number of threads to use. 
default: 4

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

1.   The reference sequence file should be in fasta format with sequence IDs containing numbers or X, Y. e.g.

>1 CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02
CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCAT
GAATCCCTAAATACCTAATTCCCTAAACCCGAAACCGGTTTCTCTGGTTGAAAATCATTGTGTATATAATGATAATTTT
ATCGTTTTTATGTAATTGCTTATTGTTGTGTGTAGATTTTTTAAAAATATCATTTGAGGTCAATACAAATCCTATTTCT
TGTGGTTTTCTTTCCTTCACTTAGCTATGGATGGTTTATCTTCATTTGTTATATTGGATACAAGCTTTGCTACGATCTA
CATTTGGGAATGTGAGTCTCTTATTGTAACCTTAGGGTTGGTTTATCTCAAGAATCTTATTAATTGTTTGGACTGTTTA
TGTTTGGACATTTATTGTCATTCTTACTCCTTTGTGGAAATGTTTGTTCTATCAATTTATCTTTTGTGGGAAAATTATT
TAGTTGTAGGGATGAAGTCTTTCTTCGTTGTTGTTACGCTTGTCATCTCATCTCTCAATGATATGGGATGGTCCTTTAG
CATTTATTCTGAAGTTCTTCTGCTTGATGATTTTATCCTTAGCCAAAAGGATTGGTGGTTTGAAGACACATCATATCAA
AAAAGCTATCGCCTCGACGATGCTCTATTTCTATCCTTGTAGCACACATTTTGGCACTCAAAAAAGTATTTTTAGATGT
TTGTTTTGCTTCTTTGAAGTAGTTTCTCTTTGCAAAATTCCTCTTTTTTTAGAGTGATTTGGATGATTCAAGACTTCTC
.........

2.   The GFF file contains 9 columns with tab-delimited. The ID for each CDS should be unique
e.g.
Chr1	TAIR9	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR9	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR9	CDS	3996	4276	.	+	2	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR9	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR9	CDS	4486	4605	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR9	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR9	CDS	4706	5095	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR9	exon	5174	5326	.	+	.	Parent=AT1G01010.1
Chr1	TAIR9	CDS	5174	5326	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR9	exon	5439	5899	.	+	.	Parent=AT1G01010.1
.......

3.   We accept two format for SNP files

I) .snp file: tab-delimited with first four columns for chr#, position, ref, and alt. e.g. the ones used for 1001 genome project 

10	2090	T	N							
10	2943	G	G						
10	2956	T	N						
10	2970	C	N						
10	3323	C	C						
10	3362	C	C						
10	4146	A	A									
10	4379	A	N

.......

II) .VCF for standard VCF file with all individuals and SNPs in one and only one file 

 
=head2 Output
Output format: Chr#,  transcript ID, total length of CDS, #of 0-folded sites, #of 4-folded sites, Pi of 0-folded sites, and Pi of 4-folded sites

Chr	Transcript	Length	L_N	L_S	Pi_N	Pi_S
1	AT1G68260.1	573	373	90	0	0	
1	AT1G06450.1	1083	706	163	0	0	
1	AT1G20560.1	1671	1092	264	1.33	1.33	
1	AT1G73170.2	1617	1018	246	1.33	0.6667	
1	AT1G17860.1	591	382	84	8	1.33	
1	AT1G75770.1	693	436	95	0	0	
1	AT1G24996.1	450	296	67	0	0	

........


=head1 Technical details



=head2 Test statistic

We calculated nucleotide diversity (Pi) over all 0-folded and 4-folded sites as well as the total number of these two kinds of sites
1) Degenerate fold for each site of the coden is calculated for each of CDS sequences in individuals sampled 
2) If different degenerate fold found between individuals, the site would not be counted in final results
3) All amibiguous sites with IUPAC DNA code will be deciphered with all possible codon with equal probability. 
   e.g. GAR (Glu) can be GAG or GAA. Thus both codon give 0,0,2 for degenerate at position 1, 2, and 3
   If codons give conflicting degenerate, follow rule 2
4) Insertion or deletion [Gap] relative to the reference is not counted
5) Sites with more than 1 alternative alleles are not counted 
6) Stop codons are not counted

=head1 AUTHORS

Jun Chen

=cut







