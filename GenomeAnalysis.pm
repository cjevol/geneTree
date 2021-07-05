package GenomeAnalysis;
use strict;
use threads;
use threads::shared;
use Thread::Queue;
use Bio::SeqIO;
use Bio::Seq;
use List::Util qw(shuffle);
use File::Copy;
use File::Basename;
use File::Find;
use Cwd qw{abs_path};
use constant GT_COL => 9; ## NOTE CHANGE TO 8 FOR STANDARD VCF....
use constant NUMH => 20; # number of header of VCF

sub gff_parser{
	print STDERR "Reading CDS info from GFF...";
	my ($gff_file, $chr_sel, $feature_ref)=@_;
	open(GFF, $gff_file) or die "Cannot open $gff_file: $!";
	my @contents=<GFF>;
	close(GFF);
	foreach (@contents){
		next if(/^#/);
		chomp;
		my ($chr, $feature, $start, $end, $strand, $attr)=(split/\t/)[0, 2,3,4,6,8];
		$chr = uc($chr);
		next unless($feature eq 'CDS');
		# if($chr=~/(\d+)/ or $chr=~/([XY])/i){ #only chromosomes with numbers and x, y
		# 	$chr=$1;
		# }else{
		# 	next;
		# }
		next unless(!$chr_sel || $chr eq $chr_sel);
		## info col can be separated by ';' or ','
		# my $parent_attr=(split/;/, $attr)[2];
		my ($parent_attr)=$attr=~/Parent=([\w|\d|\.]+)/;
		die "undefined parent id $attr" unless(defined $parent_attr);
		
		# $parent_attr=~s/Parent=mRNA\://;
		my ($transcript_id)=(split/,/, $parent_attr)[0];
		## sometimes bad annotation happens when a transcript annotated in both strands. remove later
		push @{$$feature_ref{$chr}->{$transcript_id}{$strand}}, [$start, $end];
		## 
	}
	print STDERR "\tDone!\n";
}

sub extract_CDS_by_coord{
	my ($ref_seqObj, $trans_coord_href, $transcript_id)=@_;
	my @strands=keys %$trans_coord_href; # some transcripts with two strands need to be removed
	if(scalar @strands > 1){
		print STDERR "$transcript_id has more than 1 strand, and is IGNORED!";
		return(undef,undef);
	}else{
		my $cds_seq=undef;
		my @base_genome_coord=();
		# join each exon in the CDS			
		my @exon_coords=@{$$trans_coord_href{$strands[0]}}; 
		# @exon_coords=reverse @exon_coords if($strands[0] eq '-'); # tomato or rice's gff exon 1 comes last .... comment this if not true
		foreach my $exon_coord (@exon_coords){
			my ($start, $end)=@$exon_coord;
			my $exon_seq=$ref_seqObj->subseq($start, $end);
			my @exon_base_pos=$start..$end;
			# if on the other strand, reverse complementary the sequences and coordinates
			if($strands[0] eq '-'){
				$exon_seq=reverse($exon_seq); # reverse
				$exon_seq=~tr/ATGCatgc/TACGTACG/; # complementary
				@exon_base_pos=reverse(@exon_base_pos);
			}
			$cds_seq.=$exon_seq;
			push @base_genome_coord, @exon_base_pos;
		}
		
		my $cds_seqObj=Bio::Seq->new(-id => $transcript_id, -seq => $cds_seq,-alphabet=>'dna');
		if(is_good_CDS($cds_seqObj)){
			return($cds_seqObj, \@base_genome_coord);
		}else{
			return(undef,undef); 
		}
	}
}

sub get_indSeq_by_refSeq_polym{
	my %ambig=(
    "R" => ["A", "G"],
    "Y" => ["C", "T"],
    "S" => ["G", "C"],
    "W" => ["A", "T"],
    "K" => ["G", "T"],
    "M" => ["A", "C"],
    "B" => ["C", "G", "T"],
    "D" => ["A", "G", "T"],
    "H" => ["A", "C", "T"],
    "V" => ["A", "C", "G"],
    "N" => ["A", "C", "G", "T"],
	'-' => ['-'],
	'?' => ['?'],
	);
	
	my ($cds_seq, $genome_coords_arr, $snp_list, $polymorphism_list, $ref_list, $ploid)=@_;
	my @ind_seq=split//,$cds_seq;
	foreach my $pos (0..$#ind_seq){
		# my $genome_cooord=$$genome_coords_arr[$_];
		if(exists $$snp_list{$$genome_coords_arr[$pos]}){
			my $cds_ref_allele=$ind_seq[$pos];
			my $snp_ref_allele=$$snp_list{$$genome_coords_arr[$pos]}->[0];
			my $strand=undef;
			my %complementary=('A' => 'T', 'G' => 'C', 'T' => 'A', 'C' => 'G', '-' => '-', '?' => '?',
							   'R' => 'Y', 'Y' => 'R', 'S' => 'S', 'W' => 'W', 'K'=>'M','M' => 'K',
							   'B' => 'V', 'D' => 'H', 'H' => 'D', 'V' => 'B', 'N' => 'N',
							   );
			# determine strand
			if($cds_ref_allele eq $snp_ref_allele){
				$strand=1;
			}elsif($cds_ref_allele eq $complementary{$snp_ref_allele}){
				$strand=-1;
			}else{
				print STDERR "Inconsistent genome and snp ref allele: cannot determine the strand\n";
				next;
			}
			
			$ind_seq[$pos]=$$snp_list{$$genome_coords_arr[$pos]}->[1]; # update sequence
			if($strand==-1){
				$ind_seq[$pos]=$complementary{$ind_seq[$pos]};
			}
			########################
			# add ref allele
			if(exists $$ref_list{$pos}){
				die "Inconsistent ref allele between individuals" unless ($$ref_list{$pos} eq $cds_ref_allele);
			}else{
				$$ref_list{$pos}=$cds_ref_allele;
			}
			
			## count alt alleles
			
			if($ind_seq[$pos]=~/[ATGC]/){
				if(exists $$polymorphism_list{$pos}->{$ind_seq[$pos]}){ 
					$$polymorphism_list{$pos}->{$ind_seq[$pos]}+=$ploid;
				}else{
					$$polymorphism_list{$pos}->{$ind_seq[$pos]}=$ploid;
				}
			}elsif(exists $ambig{$ind_seq[$pos]}){ 
				my $alleles=$ambig{$ind_seq[$pos]};
				# this would be wrong if not diploid, and if more than 1 alternative alleles found, the site will not be counted in later step 
				foreach my $allele (@$alleles){
					if($allele=~/\-/ || $allele=~/\?/){ # we don't count insert or deletion
						next;
					}else{ 
						if(exists $$polymorphism_list{$pos}->{$allele}){ 
							$$polymorphism_list{$pos}->{$allele}++;
						}else{
							$$polymorphism_list{$pos}->{$allele}=1;
						}
					}
				}
			}else{
				print STDERR "Found non-exist IUPAC nucleotide $ind_seq[$pos] \n"; # make sure no missing
			}
		}
	}
	return(join('', @ind_seq));
}

sub get_degenerate_fold_by_seq{
	my %nonambig=(
	 'A' => [0, 0, 4],
	 'G' => [0, 0, 4],
	 'P' => [0, 0, 4],
	 'T' => [0, 0, 4],
	 'V' => [0, 0, 4],
	 'I' => [0, 0, 3],
	 'N' => [0, 0, 2],
	 'D' => [0, 0, 2],
	 'C' => [0, 0, 2],
	 'Q' => [0, 0, 2],
	 'E' => [0, 0, 2],
	 'H' => [0, 0, 2],
	 'K' => [0, 0, 2],
	 'F' => [0, 0, 2],
	 'Y' => [0, 0, 2],
	 'W' => [0, 0, 0],
	 'M' => [0, 0, 0], 
	 '*' => [-1,-1,-1], # stop 
	);
	my %ambig=(
	 'CGT' => [0, 0, 4], # Arg
	 'CGC' => [0, 0, 4],
	 'CGA' => [2, 0, 4],
	 'CGG' => [2, 0, 4],
	 'AGA' => [2, 0, 2],
	 'AGG' => [2, 0, 2],
	 'TTA' => [2, 0, 2], # Leu
	 'TTG' => [2, 0, 2],
	 'CTT' => [0, 0, 4],
	 'CTC' => [0, 0, 4],
	 'CTA' => [2, 0, 4],
	 'CTG' => [2, 0, 4],
	 'TCT' => [0, 0, 4], # Ser
	 'TCC' => [0, 0, 4],
	 'TCA' => [0, 0, 4],
	 'TCG' => [0, 0, 4],
	 'AGT' => [0, 0, 2],
	 'AGC' => [0, 0, 2],
	);
	
	my ($seqObj)=@_;
	my $dna_seq=uc($seqObj->seq);
	my $aa_seq=$seqObj->translate->seq;
	my $aa_lgth=length($aa_seq);

	my @degeneracy_fold=();
	for(my $i=0; $i<$aa_lgth;$i++){
		my $aa=substr($aa_seq, $i,1);
		if(exists $nonambig{$aa}){
			push @degeneracy_fold, @{$nonambig{$aa}};
		}else{
			my $genetic_codon=substr($dna_seq, $i*3, 3);
			if(exists $ambig{$genetic_codon}){ # aa eq 'R', 'L' or 'S' 
				push @degeneracy_fold,  @{$ambig{$genetic_codon}};
			}elsif($genetic_codon=~/\-/ || $genetic_codon=~/\?/){ ## contains insertion or deletion
				push @degeneracy_fold, (-1,-1,-1);
			}elsif($genetic_codon=~/[^ATGC]/){ # must contain IUPAC code eg, "K","R","W",'N' ....et ac 
				my @codon_combinations=get_allFree_combinations($genetic_codon);
				my @codon_fold=();
				my @final_codon_fold=();
				foreach (@codon_combinations){
					my @fold=get_degenerate_fold_by_seq(Bio::Seq->new(-seq=>$_,-alphabet=>'dna'));
					$codon_fold[0]{$fold[0]}=1;
					$codon_fold[1]{$fold[1]}=1;
					$codon_fold[2]{$fold[2]}=1;
				}
				foreach (@codon_fold){
					if(scalar keys %$_ == 1){ # for each position only one degeneracy fold shall be assign, skip if otherwise
						push @final_codon_fold,keys %$_;
					}else{
						push @final_codon_fold,-1;
					}
				}
				push @degeneracy_fold, @final_codon_fold;	
			}else{
				push @degeneracy_fold, (-1,-1,-1);
				print STDERR "\t\tUnkown codon found $genetic_codon!\n"; # That error shall not happen but still check
			}
		}
	}
	return(@degeneracy_fold);
}
sub get_allFree_combinations{
	my ($genetic_codon)=@_;
	my @bases=split//,$genetic_codon;
	my @all_combinations=();
	my %ambig=(
    "R" => ["A", "G"],
    "Y" => ["C", "T"],
    "S" => ["G", "C"],
    "W" => ["A", "T"],
    "K" => ["G", "T"],
    "M" => ["A", "C"],
    "B" => ["C", "G", "T"],
    "D" => ["A", "G", "T"],
    "H" => ["A", "C", "T"],
    "V" => ["A", "C", "G"],
    "N" => ["A", "C", "G", "T"],
	'A' => ['A'],
	'T' => ['T'],
	'G' => ['G'],
	'C' => ['C'],
	);
	my $first_base_arr=$ambig{$bases[0]};
	my $second_base_arr=$ambig{$bases[1]};
	my $third_base_arr=$ambig{$bases[2]};
	foreach my $base1 (@$first_base_arr){
		foreach my $base2 (@$second_base_arr){
			foreach my $base3 (@$third_base_arr){
				push @all_combinations, $base1.$base2.$base3;
			}
		}
	}	
	return(@all_combinations);
}


# sub parallele_machine_bundle{
# 	my ($array, $num_threads, $function_ref, @options)=@_;
# 	my $num_total_jobs=scalar @$array;
# 	my $num_jobs_per_child=int $num_total_jobs/$num_threads;
# 	die "Distribute $num_jobs_per_child transcripts among $num_threads threads?" unless($num_jobs_per_child); # we assume more transcripts than threads ...
# 	my @child_process=();
# 	push @child_process, [splice @$array, 0, $num_jobs_per_child] for(1..$num_threads);
# 	push @$child_process[-1],@$array;
# 	parallele_machine(\@child_process, $function_ref, @options);
# }
#

sub process_array_ithread{
	my ($arr_ref, $num_threads, $func_ref, @options)=@_;
	my @threads;
	my $q=Thread::Queue->new();
	push @threads, async{$func_ref->($q, @options)} for(1..$num_threads);	
	$q->enqueue($_) for @$arr_ref;
	$q->enqueue(undef) for @threads;
	$_->join() for @threads;
}

sub parallele_machine{
	my ($array, $function, @options)=@_;
	my @child;
	foreach my $child_process (@$array){
		my $pid=fork();
		if($pid){
			push @child, $pid; # parent process
		}elsif($pid==0){ # child process
			$function->($child_process, @options);
			exit 0;
		}else{
			die "Cannot fork: $!";
		}
	}
	foreach (@child){
		my $tmp=waitpid($_,0);		
	}
}

sub parallele_machine_by_splice{
	my ($array, $num_threads, $function_ref, @options)=@_;
	my $array_length=scalar @$array;
	my @child_process=();
	# initiate
	if($num_threads <= $array_length){
		@child_process=splice(@$array, 0, $num_threads);
	}else{
		@child_process=@$array;
	}
	
	while(@child_process){
		parallele_machine(\@child_process, $function_ref, @options);
		if($num_threads <= $array_length){
			@child_process=splice(@$array, 0, $num_threads);
		}else{
			@child_process=splice(@$array, 0);
		}
	}
}

sub split_snp_file{
	print STDERR "Sample and split snp files by chromosome\n";
	my ($snpdir, $snp_type, $sample_size, $snp_tmp_dir, $chr_ids_ref, $num_threads)=@_;
	# source_dir, $type, $sample_size, $target_dir, $chr, $num_proc
	my @sampled_filenames=sample_files_by_size($snpdir, $snp_type, $sample_size);
	split_snp_by_chromome_parallele(\@sampled_filenames, $snp_tmp_dir, $chr_ids_ref, $num_threads);
	my @snp_filenames=get_filenames_recursive($snp_tmp_dir, $snp_type);
	print STDERR "\tDone!\n";
	return(@snp_filenames);
}

sub split_vcf_file{
	#####################
	my ($snpdir, $sample_size, $snp_tmp_dir, $chr_ids_ref, $num_threads)=@_;
	my @sample_filenames=get_filenames_recursive($snpdir, 'vcf');
	die "We only accept ONE file in VCF mode now" unless (scalar @sample_filenames==1);
	my $ind_name=fileparse($sample_filenames[0]);
	print STDERR "Sample and split the VCF file $ind_name by chromosome\n";
	
	### assign filehandles to array
	my %fh=();
	foreach my $key (@$chr_ids_ref){
		my $target_filename=$snp_tmp_dir.'/'.$key.'-'.$ind_name;
		open(my $fh, '>', $target_filename) or die "Cannot write into $target_filename: $!";
		$fh{$key}=[$fh,0];
	}	
	open(VCF, $sample_filenames[0]) || die "Cannot open $sample_filenames[0]: $!";
	my @info_lines = map~~<VCF>, 1..NUMH;
	my $header_line=<VCF>;
	chomp $header_line;
	my @header_cols=split/\t/,$header_line; 
	my @ind_ids=@header_cols[GT_COL..$#header_cols]; #get the ids
	my @ind_indice=random_permuted_index(scalar @ind_ids);
	@ind_indice=@ind_indice[0..($sample_size-1)];
	print STDERR "$sample_size samples have been selected by random: ",join("\t", @ind_ids[@ind_indice]),"\n";
	
	$_+=GT_COL for @ind_indice;
	my @selected=0..(GT_COL-1);
	push @selected, @ind_indice;
	
 	while(<VCF>){
		chomp;
		my @cols=split/\t/;
		my $chr=$cols[0];
		# if($chr=~/(\d+)/ or $chr=~/([XY])/i){
		# 	$chr=$1;
			if(exists $fh{$chr}){
				if($fh{$chr}[1]==0){
					$fh{$chr}[1]=1;
					print {$fh{$chr}[0]} @info_lines;
					print {$fh{$chr}[0]} join("\t", @header_cols[@selected]),"\n"; 
				}
				print {$fh{$chr}[0]} join("\t", @cols[@selected]),"\n";
			}# else{
# 				die "Cannot open proper filehandle to write for chromosome $chr: $!";
# 			}
		# }
	}
	close(VCF);
	## close filehandles
	foreach (keys %fh){
		close $fh{$_}[0];
	}
	print STDERR "\tDone!\n";
	return(@sample_filenames);
}

sub split_snp_by_chromome_parallele{
	print STDERR "Splitting SNP files by chromomsome...";
	my ($source_filenames_aref, $snp_tmp_dir, $chr_ids_aref, $num_threads)=@_;
	parallele_machine_by_splice($source_filenames_aref, $num_threads, \&split_file_by_1stCol, ($snp_tmp_dir, $chr_ids_aref));	
}


sub split_file_by_1stCol{
	####################################
	#  split into each chromosome
	####################################
	
	my ($source_file, $target_dir, $keyword_aref)=@_;
	my $ind_name=fileparse($source_file);
	my %fh=();
	### assign filehandles to array
	foreach my $key (@$keyword_aref){
		my $target_filename=$target_dir.'/'.$key.'-'.$ind_name;
		open(my $fh, '>', $target_filename) or die "Cannot write into $target_filename: $!";
		$fh{$key}=$fh;
	}
	## split by first column chromosome
	open(IN, $source_file) or die "Cannot open $source_file: $!";
	while(<IN>){
		my ($chr)=(split/\t/)[0];
		$chr=uc($chr);
		# if($chr=~/(\d+)/ or $chr=~/([XY])/i){
		# 	$chr=$1;
			if(exists $fh{$chr}){
				print {$fh{$chr}} $_;
			}# else{
# 				die "Cannot open proper filehandle to write for chromosome $chr: $!";
# 			}
		# }
	}
	close(IN);
	## close filehandles
	foreach (keys %fh){
		close $fh{$_};
	}
}

sub sample_files_by_size{
	########################################
	#  sample x files from source directory
	########################################
	my ($source_dir, $file_type, $sample_size)=@_;
	my @source_filenames=get_filenames_recursive($source_dir, $file_type);
	die "No $file_type file exists in $source_dir" unless(@source_filenames);	
	my @sampled_filenames=sample_array($sample_size, @source_filenames);
	print STDERR "We have been sampled $sample_size files: \n\t", join("\n\t", @sampled_filenames), "\n";
	return(@sampled_filenames);
}

sub sample_array{
	# pick n samples from a list
	my ($size, @array)=@_;
	my $array_size=scalar @array;
	if($array_size<$size){
		$size=$array_size;
		print STDERR "Warning: reduce the sample size $size to $array_size!\n"; 
	}
	$size--;
	my @premuted_array=random_permutation(@array);
	return(@premuted_array[0..$size]);
}

sub random_permutation{
	# shuffle the array
	my (@array)=@_;
	my $size=scalar @array;
	my @permuted_index=random_permuted_index($size);
	return(@array[@permuted_index]);
}

sub random_permuted_index{
	# return permuted index in a array
	my ($n)=@_;
	$n--;	
	return(shuffle(0..$n));
}

sub is_good_CDS{
	my ($cds_seqObj)=@_;
	my $id=$cds_seqObj->id;
	if($cds_seqObj->length % 3!=0){
		print STDERR "\t\t$id is not properly annotated, and is IGNORED!\n";
		return(0);
	}else{
		my $aa_seq=$cds_seqObj->translate->seq;
		if($aa_seq!~/^M.+\*$/){
			print STDERR "\t\t$id does not start with 'ATG',or end with stop!\n";
			return(1); #choose the action with partial cds
		}elsif($aa_seq=~/\w+\*\w+/){
			print STDERR "\t\t$id contains premature stop codon, and is IGNORED!\n";
			return(0);
		}else{
			return(1);
		}
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

sub get_filenames_recursive{
	my ($parent_path, $filter)=@_;
	$parent_path=abs_path($parent_path);
	my @filenames=();
	find sub{
		push @filenames, $File::Find::name if(-f $File::Find::name && grep(/\.$filter$/i, $_));
	}, $parent_path;
	return(@filenames);
}

1