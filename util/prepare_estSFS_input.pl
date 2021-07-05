#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
use File::Copy;
use File::Basename;
use File::Find;
use Cwd qw{abs_path};


my $ingr_path = shift;
my $outgr_path = shift;
my %outgr=();
my $outgr_index=0;
foreach (@ARGV){
	$outgr_index++;
	my @inds = split/,/;
	$outgr{$_}=$outgr_index foreach (@inds);
}
foreach (sort {$outgr{$a} <=> $outgr{$b}} keys %outgr){
	print STDERR $_,"\t", $outgr{$_},"\n";
}

my %ingr_freq =();
print STDERR "Reading focus alignments $ingr_path\n";
get_freq_aln($ingr_path, \%ingr_freq);

print STDERR "\nReading outgr alignments $outgr_path and printing\n";
foreach my $seq_name (keys %ingr_freq){
	my %outgr_freq =();
	my @snp_pos = sort {$a <=> $b} keys %{$ingr_freq{$seq_name}};
	get_outgr_states($outgr_path, $seq_name, \@snp_pos, \%outgr_freq, \%outgr);
	foreach my $pos (@snp_pos){
		print "$seq_name\t$pos\t";
		if(exists $ingr_freq{$seq_name}->{$pos}{'A'}){
			print $ingr_freq{$seq_name}->{$pos}{'A'},",";
		}else{
			print "0,";
		}
		if(exists $ingr_freq{$seq_name}->{$pos}{'C'}){
			print $ingr_freq{$seq_name}->{$pos}{'C'},",";
		}else{
			print "0,";
		}
		if(exists $ingr_freq{$seq_name}->{$pos}{'G'}){
			print $ingr_freq{$seq_name}->{$pos}{'G'},",";
		}else{
			print "0,";
		}
		if(exists $ingr_freq{$seq_name}->{$pos}{'T'}){
			print $ingr_freq{$seq_name}->{$pos}{'T'};
		}else{
			print "0";
		}
		
		foreach my $gr (1..$outgr_index){
			if(exists $outgr_freq{$gr}){
				if(exists $outgr_freq{$gr}->{$pos}{'A'}){
					print "\t1";
				}else{
					print "\t0";
				}
			
				if(exists $outgr_freq{$gr}->{$pos}{'C'}){
					print ",1";
				}else{
					print ",0";
				}

				if(exists $outgr_freq{$gr}->{$pos}{'G'}){
					print ",1";
				}else{
					print ",0";
				}
			
				if(exists $outgr_freq{$gr}->{$pos}{'T'}){
					print ",1";
				}else{
					print ",0";
				}
			}else{
				print "\t0,0,0,0";
			}
		}
		print "\n";
	}
}
print STDERR "Done\n";

#########################
sub get_outgr_states{
	my %parse_gt=(
		'A' => 'AA',
		'T' => 'TT',
		'G' => 'GG',
		'C' => 'CC',
		'R' => 'AG',
		'Y' => 'CT',
		'S' => 'GC',
		'W' => 'AT',
		'K' => 'GT',
		'M' => 'AC'
	);
	
	my ($path, $seq_name, $snp_pos, $hash, $gr) = @_;
	my $file = $path.'/'.$seq_name.'.fas';
	my $in = Bio::AlignIO->new(-file => $file, -format => 'fasta');
	my $aln = $in->next_aln;
	foreach my $snp_pos (@$snp_pos){
		my $aln2 = $aln->slice($snp_pos, $snp_pos);
		foreach my $seqObj ($aln2->each_seq){
			my $id = $seqObj->id;
			my $group_id = undef;
			if(exists $$gr{$id}){
				$group_id = $$gr{$id};
			}else{
				next;
			}
			
			my $gt = $seqObj->seq;
			next if($gt eq 'N');
			
			my @alleles=();
			if(exists $parse_gt{$gt}){
				@alleles = split//,$parse_gt{$gt};
			}else{
				die "Unknown gt $seq_name, $snp_pos, $gt ";
			}
			
			foreach my $base (@alleles){
				if(exists $$hash{$group_id}->{$snp_pos}{$base}){
					$$hash{$group_id}->{$snp_pos}{$base}++;
				}else{
					$$hash{$group_id}->{$snp_pos}{$base}=1;
				}	
			}
			my @base = sort {$$hash{$group_id}->{$snp_pos}{$b} <=> $$hash{$group_id}->{$snp_pos}{$a}} keys %{$$hash{$group_id}->{$snp_pos}};
			if(scalar @base > 1){
				shift @base;
				delete $$hash{$group_id}->{$snp_pos}{$_} foreach(@base);
			}
		}	
	}
}



sub get_freq_aln{
	my %parse_gt=(
		'A' => 'AA',
		'T' => 'TT',
		'G' => 'GG',
		'C' => 'CC',
		'R' => 'AG',
		'Y' => 'CT',
		'S' => 'GC',
		'W' => 'AT',
		'K' => 'GT',
		'M' => 'AC'
	);
	my ($path, $hash) = @_;
	my @files = get_filenames_recursive($path, 'fas');
	my $countor = 0;
	foreach my $file (@files){
		$countor++;
		print STDERR "." if($countor % 1000 == 0);
		my $seq_name = basename($file, '.fas');
		my $in = Bio::AlignIO->new(-file => $file, -format => 'fasta');
		my $aln = $in->next_aln;
		my $iupac = $aln->consensus_iupac();
		# print "$seq_name\n",$iupac,"\n";
		my @iupac = split//, $iupac;
		my @snp_pos = ();
		foreach (0...$#iupac){
			push @snp_pos, $_ if($iupac[$_] ne 'A' | $iupac[$_] ne 'T' | $iupac[$_] ne 'G' | $iupac[$_] ne 'C' );
		}
		
		foreach my $snp_index (@snp_pos){
			my $snp_pos = $snp_index + 1;
			my $aln2 = $aln->slice($snp_pos, $snp_pos);
			foreach my $seqObj ($aln2->each_seq){
				my $gt = $seqObj->seq;
				next if($gt eq 'N');
				my @alleles=();
				if(exists $parse_gt{$gt}){
					@alleles = split//,$parse_gt{$gt};
				}else{
					die "Unknown gt $seq_name, $snp_pos, $gt ";
				}
				foreach my $base (@alleles){
					if(exists $$hash{$seq_name}->{$snp_pos}{$base}){
						$$hash{$seq_name}->{$snp_pos}{$base}++;
					}else{
						$$hash{$seq_name}->{$snp_pos}{$base}=1;
					}	
				}
			}
			if(scalar keys %{$$hash{$seq_name}->{$snp_pos}} ==1){
				delete $$hash{$seq_name}->{$snp_pos};
			}
		}
	}
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
