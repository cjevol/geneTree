#!/usr/bin/perl -w
use strict;
use File::Copy;
use File::Basename;
use File::Find;
use Cwd qw{abs_path};

my $configure = shift;
print STDERR "##### 1 configure file ############################\n";
open(IN, $configure) or die "Lack configure file: $!";
my %options=();
while(<IN>){
	chomp;
	next if(/^\s*\#*$/ or /^#/);
	my ($key, $value) = (split/\=/)[0,1];
	$key=~s/\s//g;
	$value=~s/\s//g;
	$value=~s/#.*$//;
	$value=~s/\'//g;
	$value=~s/\"//g;
	print "$key\t$value\n";
	$options{$key}=$value;
}
close(IN);
## get genome fasta
my ($ref_file) = get_filenames_recursive($options{'genome_dir'}, qw/fa.gz fas.gz fasta.gz fa fas fasta/);
my ($gff_file) = get_filenames_recursive($options{'genome_dir'}, qw/gff gff3 gtf gff.gz gff3.gz gtf.gz/);
if(-e $ref_file){
	print STDERR "found reference genome $ref_file\n";
}else{
	die "Lack ref fasta file";
}

if(-e $gff_file){
	print STDERR "found gff file $gff_file\n";
}else{
	die "Lack gff file";
}

# step convert vcf2snp
{
	print STDERR "\n##### 2 convert vcf file ############################\n";

	{ ## comment this block to repeat without change outgr 
		mkdir 'tmp_snp' or die "Cannot create tmp_snp folder" unless(-e 'tmp_snp');
		my $cmd = "perl util\/convert-vcf2snp.pl $gff_file $options{'vcf_dir'} tmp_snp ";
		print STDERR "$cmd\n";
		system($cmd) == 0 or die "failed to convert snp: $!";
	}

	mkdir 'snp' or die "Cannot create snp folder" unless(-e 'snp');
	if($options{'Pop_tags'}=~/,/){
		my @pop_files = split/,/,$options{'Pop_tags'};
		foreach (@pop_files){
			my $file = 'tmp_snp/'.$_.'.snp';
			my $cmd = "cp $file snp";
			system($cmd) == 0 or die "failed to mv $_: $!";
		}
	}else{
		my $cmd = "cp tmp_snp/".$options{'Pop_tags'}.'*.snp snp';
		system($cmd) == 0 or die "failed to mv snp files: $!";
	}
}

# ## calculate piN/piS rate
print STDERR "\n##### 3 Getting polymorphic sites ############################\n";

my $output_file = 'piNpiS_Ind'.$options{'num_inds'};
$output_file.='_diploid' if($options{'ploid'} == 2);
$output_file.='.txt';
my $log_file = 'log_'.$output_file;
#
{
	my $cmd = "perl calc-piNpiS-ithreads.pl -r $ref_file -g $gff_file -p snp -s $options{'num_inds'} -d $options{'ploid'} -t $options{'threads'} -o 1 -w focus_tmp 1> $output_file 2> $log_file";
	print STDERR "$cmd\n";
	system($cmd) == 0 or die "fail to run calc-piNpiS-ithreads.pl: $!";
}

#
# ##

if($options{'unfolded'}){
	print STDERR "\n##### 4 Polarize frequency ############################\n";
	my $num_outgr=0;
	{
		print STDERR "Get outgr aln\n";

		if(-d 'outgr_snp'){
			# system('rm -r outgr_snp') or "fail to delete outgr_snp: $!";
		}else{
			system("mkdir outgr_snp") ==0 or die "fail to create folder outgr_snp: $!";
		}
		if(exists $options{'Outgr1'}){
			my @files = split/,/,$options{'Outgr1'};
			foreach (@files){
				my $outgr1='tmp_snp/'.$_.'.snp';
				system("cp $outgr1 outgr_snp") ==0 or die "fail to mv: $!";
			}
			$num_outgr++;
		}

		if(exists $options{'Outgr2'}){
			my @files = split/,/,$options{'Outgr2'};
			foreach (@files){
				my $outgr2='tmp_snp/'.$_.'.snp';
				system("cp $outgr2 outgr_snp") ==0 or die "fail to mv: $!";
			}
			$num_outgr++;
		}
		if(exists $options{'Outgr3'}){
			my @files = split/,/,$options{'Outgr3'};
			foreach (@files){
				my $outgr3='tmp_snp/'.$_.'.snp';
				system("cp $outgr3 outgr_snp") ==0 or die "fail to mv: $!";
			}
			$num_outgr++;
		}
		
		my $cmd = "perl util\/get_aln_from_snp.pl -r $ref_file -g $gff_file -p outgr_snp -s $num_outgr -d $options{'ploid'} -t $options{'threads'} -o 1 -w outgr_tmp 1>/dev/null 2>&1 ";
		print STDERR $cmd,"\n";
		system($cmd) == 0 or die "fail to get outgr aln: $!";	
	}

	{
		print STDERR "# prepare input (take a while) \n";
		my $outgr_str = $options{'Outgr1'};
		$outgr_str.=' '.$options{'Outgr2'} if(exists $options{'Outgr2'});
		$outgr_str.=' '.$options{'Outgr3'} if(exists $options{'Outgr3'});
		my $cmd = "perl util\/prepare_estSFS_input.pl focus_tmp\/align\/ outgr_tmp\/align\/ $outgr_str > estSFS_full.txt";
		system($cmd) == 0 or die "fail to generate estSFS input file: $!";
	}
	
	## remove some sites not in frequency.txt to save time
	{
		my $cmd = "perl util\/filter_estSFS.pl estSFS_full.txt frequency.txt";
		system($cmd) == 0 or die "fail to filter estSFS input file: $!";
	}
	
	{
		open(IN, "estSFS_full.txt") or die "Cannot open estSFS_full.txt: $!";
		open(OUT,">estSFS_slim.txt");
		my @content=<IN>;
		chomp @content;
		foreach (@content){
			my @cols = split/\t/;
			print OUT "$cols[2]";
			print OUT "\t$cols[$_+2]" for(1..$num_outgr);
			print OUT "\n";
		}
		close(IN);
		close(OUT);
	}

	if($num_outgr == 1){
		my $cmd = "perl util\/find_fixation.pl estSFS_full.txt > anc_info.txt";
		system($cmd) == 0 or die "fail to parse anc result: $!";
	}else{
		print STDERR "# running estSFS \n";
		open(CON, ">config-rate6.txt");
		print CON "n_outgroup $num_outgr\nmodel 2\nnrandom 10\n";
		close(CON);
		my $cmd = "util\/est-sfs config-rate6.txt estSFS_slim.txt util\/seedfile.txt out_sfs out_p_anc ";
		system($cmd) == 0 or die "fail to run estSFS: $!";
		# parse results
		$cmd = "perl util\/parse_estSFS_output.pl out_p_anc estSFS_full.txt $num_outgr > anc_info.txt";
		system($cmd) == 0 or die "fail to parse estSFS result: $!";
	}

	{
		my $cmd = "perl util\/polarize_freq.pl frequency.txt anc_info.txt > frequency_DAF.txt";
		system($cmd) == 0 or die "fail to polarize frequency.txt: $!";
		$cmd = "perl util\/parse_freq.pl frequency_DAF.txt";
		system($cmd) == 0 or die "fail to parse frequency: $!";
	}

}else{
	my $cmd = "perl util\/parse_freq.pl frequency.txt";
	system($cmd) == 0 or die "fail to parse frequency: $!";

}

system("perl util\/summaryFreqs-by-fold2.pl frequency.txt > num_snps.txt") == 0 or die "Fail to summarize freq file";
## mv results
system("mkdir $options{'out_dir'}")==0 or die "fail to build output folder: $!";
if($options{'unfolded'}){
	my $cmd = "mv out_sfs out_p_anc anc_info.txt frequency.txt frequency_DAF.txt freq_0.txt freq_4.txt missing_ratio.txt num_snps.txt $output_file $log_file $options{'out_dir'}";
	system($cmd) == 0 or die "fail to move results: $!";
}else{
	my $cmd = "mv frequency.txt freq_0.txt freq_4.txt missing_ratio.txt num_snps.txt $output_file $log_file $options{'out_dir'}";
	system($cmd) == 0 or die "fail to move results: $!";
}
#
## Run R script
{
	print STDERR "\n##### 5 Calculate piN/piS ############################\n";
	open(ROUT,">Rscript.r");
	print ROUT "source\('util\/Rfunctions.r'\)\n";
	my $Roptions= "rdir=\"$options{'genome_dir'}\", outdir=\"$options{'out_dir'}\"";
	$Roptions.= ', filter='.$options{'filter'} if(exists $options{'filter'});
	$Roptions.= ', allSite='.$options{'allSite'} if(exists $options{'allSite'});
	$Roptions.= ', common.list="'.$options{'geneList'}.'"' if(exists $options{'geneList'});
	$Roptions.= ', write.genes="'.$options{'write_genes'}.'"' if(exists $options{'write_genes'});
	$Roptions.= ', size='.$options{'size'} if(exists $options{'size'});

	print ROUT "res=calculate_piNpis\($Roptions\)\n";
	print ROUT "write.table(as.data.frame(t(res)),file='piNpiS_Rout.txt',quote=F,row.names=F,col.names=T,sep=\"\\t\")";
	my $cmd="R CMD BATCH Rscript.r";
	system($cmd) == 0 or die "fail to move results: $!";
	system("mv Rscript.r.Rout $options{'out_dir'}") == 0 or print STDERR "fail to move Rscript.r.Rout\n";

	print STDERR "Done!\n";
	print STDERR "#################################\n";
}


#######################################
sub get_filenames_recursive{
	my ($parent_path, @filter)=@_;
	$parent_path=abs_path($parent_path);
	my @filenames=();
	foreach my $filter (@filter){
		find sub{
			push @filenames, $File::Find::name if(-f $File::Find::name && grep(/\.$filter$/i, $_));
		}, $parent_path;
	}
	return(@filenames);
}
