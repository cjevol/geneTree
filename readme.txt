(https://github.com/cjevol/geneTree/tree/main)


0) The pipleline will read SNPs one chromosome each time and spreads to multiple thread to avoid memory overflow 
However, if the fasta file contains thousands of scaffolds / contigs instead of chromosomes, this will cause the system to complain for opening too many filehandles. 
One can try to increase the number of open files:
ulimit -n 100000

1) requirements 
Perl modules:
Bioperl module Bio::Seq and Bio::SeqIO (see https://github.com/bioperl/bioperl-live/blob/master/README.md for installation)
(other perl modules: threads, Thread, Time, Benchmark, etc. should be available for most linux system)
estSFS (if two or more outgroups available. to download http://www.homepages.ed.ac.uk/pkeightl/software.html)


2) input files (Please put fasta and gff files in one folder)

	i) reference genome file in fasta format (xx.fasta/ fa/ fas) 

	ii) gene annotation file in GFF format (xx.gff/gff3)
	
	iii) vcf.gz file(s) (xx.vcf.gz, can be multiple files in one folder)
	
# please note gff file should be sorted (exon 1 comes first even if on anti-sense strand)
# rows contain 'CDS' info will be extracted 
# Parent field in last column will be extracted for transcript ID (only word, dot, and numbers are allowed otherwise change the match pattern in GenomeAnalysis.pm line 36)
# e.g.: 
# Hic_asm_11      HR      CDS     32016767        32016806        .       -       0       ID=cds.evm.model.HIC_ASM_11.2140;Parent=evm.model.HIC_ASM_11.2140
# Hic_asm_11      HR      CDS     32016572        32016627        .       -       2       ID=cds.evm.model.HIC_ASM_11.2140;Parent=evm.model.HIC_ASM_11.2140
# Hic_asm_11      HR      CDS     32014987        32015039        .       -       0       ID=cds.evm.model.HIC_ASM_11.2140;Parent=evm.model.HIC_ASM_11.2140
# .....


3) configure file
A file contains all parameters (space will be ignored and lines begin with # are comments)

one can put vcf.gz files in 'vcf_dir' folder; fasta and gff files in 'genome_dir' folder

If outgroup(s) exists in vcf, one can set 'unfolded = 1' (otherwise 0) and specifies names of outgroups.
At most three outgroups are allowed with with Outgr1 the closet and Outgr3 the most distant in genetic distance (see manual of estSFS)
Each Outgr can include multiple individuals delimited by comma (Outgr1 = 'LP197014', or Outgr1 = 'LP197014,LP197015,LP197016' if they are from same species)  
To run estSFS (Keightley) one must have at least two outgroups.  With only one outgroup fixation sites of Outgr1 will be assumed as ancestral states

"Pop_tag" can be a string that pop ID starts with (e.g. Pop_tag='BWL' for BWL01, BWL02, etc) or individual IDs delimited by comma (Pop_tag='BWL01,BWL02,BWL03') 

"num_inds" is number of individuals to sample (if fewer than the actual number of individuals specified by Pop_tag in VCF, random downsample will be used)

"threads" is number of threads used to extract CDS regions and summarize frequency. (4 or 5 threads should be enough and it won't take too long for most species. Too many threads also means more RAM) 

"filter, allSite, geneList, size and write_genes" are optimal options for estimating piN/piS ratio 
"filter" (T/F. default: T) whether to perform filtering on quality of genes. (by default genes of missing rate larger than 0.5 or more than 10 missing sites will be removed)
filtering will also be performed on blastp criteria (eval<1e-20, bscore>200, query coverage > 90%, query length>30 a.a) if a file contains keyword "blastp" exists in "genome_dir" folder
	e.g. 
		 Transcript qlen subject slen eval score identity
		 evm.model.ORIGINAL_SCAFFOLD_277.4	107	AP1S1_ARATH	161	5.83e-37	122	87
		 evm.model.FRAGSCAFF_SCAFFOLD_67.17	124	MATK_BARAL	510	4.1	26.9	53
one can gain the blastp file by blast protein sequences agains Swissprot database (blastp -query xx -db xx -out blastp.txt -outfmt '6 qseqid qlen sseqid slen evalue bitscore pident'). Of course blastn with nucleotide search is similar if one prefers
the default filtering cutoffs can be changed in 'calculate_piNpis()' function of Rfunction.r file

filtering will also be performed on longest CDS. If multiple coding models (alternative splicing) exist for the same gene, one can choose to use only the longest (or whatever criteria one prefers). Just put gene names selected in one file named "longest.txt" with a header "Transcript" and put the file in "genome_dir" folder (i.e. an include list with filtering)
	e.g.	
		Transcript
		evm.model.ORIGINAL_SCAFFOLD_277.4
 		evm.model.FRAGSCAFF_SCAFFOLD_67.17
 		.....

filtering will be performed on completeness of CDS (containing both start and stop codon). If one needs to remove some genes (missing start or stop codon or containing premature stop codon), one can create a file named 'partial.txt' with a header "Transcript" and put the file in "genome_dir" folder (i.e. an exclude list)
	e.g.	
		Transcript
		evm.model.ORIGINAL_SCAFFOLD_277.4
 		evm.model.FRAGSCAFF_SCAFFOLD_67.17

Please note that if filtering is on, we will perform a resample procedure to make sure same ratio of nonpolymorphic genes to polymorphic genes. This means that the final results may have slightly fewer genes. (one can turn off this by change resample=F in 'calculate_piNpis()' function in Rfunction.r) 


"size" number of genes to downsample to

"geneList" another gene list to include with a header "Transcript" (calculation will be done only for those genes and NO FILTERING will be performed on this set of genes)
	e.g.
		Transcript
		evm.model.ORIGINAL_SCAFFOLD_277.4
 		evm.model.FRAGSCAFF_SCAFFOLD_67.17

"write_genes" if you want to see detail of piN/piS for each genes

"allSite" (T/F,default:F) if nonvariant (nonpolymorphic sites) were also included in VCF. 
Please note it is important to know how many sites were actually sequenced for each gene including polymorphic and non-polymorphic sites. But usually one only keeps SNP sites in VCF. In this case the length of non-poly sites will be estimated from the same missing ratio as in SNP sites.


###################### an example of contig.txt 

vcf_dir = "/data/cjevol/test/input/"						# vcf directory
genome_dir = '/data/cjevol/test/genome/'  					# where keeps gff and fasta
out_dir = '/data/cjevol/test/output/'						# output directory
unfolded = 1 												# 1:unfolded sfs; 0: folded
Pop_tags  = 'BWL'  											# Pop tag or Individual IDs separated by comma
Outgr1 = 'LP197014' 										# outgr1 inds separated by comma ('LP197014,LP197015,LP197016')
Outgr2 = 'LP197068' 										# # outgr2 inds separated by comma
Outgr3 = 'TMS17' 											# # outgr3 inds separated by comma
num_inds = 10 												# number of individuals to select
ploid = 2 													# ploid
threads = 4 												# number of threads 
# filter = T 													# to filter (T/F) genes with quality in R script
# allSite = F													# if vcf contains non-polymorphic sites as well
# geneList = 'filename'										# if you have a list of genes to calculate piN/piS on 
# size = 1000 												# downsample the data to a size of 1000 
# write_genes = 'output_piNpiS_details' 					# output a detail list of piNpiS for each genes


4) estSFS (download http://www.homepages.ed.ac.uk/pkeightl/software.html)
We use estSFS by keightley and Jackson (2018) to polarize SFS. The program asks for no more than 200 alleles in total and can only handle around 1 M SNP at a time. So I only polarize SNPs in CDS regions which in most case should be fewer than that. 
However, if one has too many SNPs, one can separate the input file "estSFS_slim.txt" and run estSFS and this pipeline step by step (or contact me). Sorry for inconvenience 

I used the rate-6 model for estSFS one can also choose 0= Jukes Cantor model or 1=Kimura 2-parameter model (line 160 model in main.pl) one may also change nrandom for each model (see util/config-kimura.txt and config-JC.txt)


5) output files
All output files will be in folder output
"piNpiS_Rout.txt" is the final output of summarized piN/piS ratio together with some other info. One can also choose to print details for each gene by uncomment "write_genes" in config.txt

"Rscript.r.Rout" also contains some log information when calculating piN/piS including possible error messages if it fails

"piNpiS_Indxx_diploid.txt" contains raw piN and piS without rescaling and filtering

"out_p_anc" output of estSFS. "anc_info.txt" parse estSFS results with coordinates

"dfe_input.txt" input of polyDFE2 if unfold SFS is chosen

"sfs.pdf" plot of SFS

and some other files may not be useful for everyone. 

6) to repeat pipeline on another populations
To re-run pipeline on another populations if same VCF and outgroups are used

	rename output folder ('mv output xx')
	
	remove some folder ('rm -r focus_tmp snp')
	
	comment main.pl line 44-49 (step 2)

	change Pop_tags in config.txt
	
	run "main.pl config.txt"

7) errors
Usually step 1 and step 2 should not run into errors. So if any errors occur in later steps, one can alway comment main.pl line 44-49 to skip step 2 to save some time 
i) step 3 should create a series of output, one can check if anything printed in "piNpiS_Indxx_diploid.txt"; If this file is empty one can check the log file "log_piNpiS_Indxx_diploid.txt" for error message. 
	The log file prints out the number of SNP files and number of SNPs found in each file. If there is '0' SNPs found it usually means a) that ID in fasta and GFF files do not match.b) GFF parent field in wrong format or with unidentified characters (one can modify line 36 in GenomeAnalysis.pm). 
	The log file may also complain if out of RAM occurs, like 'Thread xx has been terminated'.  
	The log file may complain "Inconsistent genome and snp ref allele: cannot determine the strand" which means site in reference fasta doesn't match ref allele in VCF. The site will be skipped. So no worries.
	The log file may complain "xx transcript is not properly annotated or does not start with 'ATG' or end with stop codon". No worry nothing will be removed. But if too many complain found one might need to check if the GFF is corrected sorted especially for those in antisense strand  
	Most of errors in this step are caused by GFF format. 
	
ii) step 4 polarize SFS
Error may be caused by too many SNPs (estSFS can deal with ~1M SNPs). One has to slit the input and run separately if more SNPs needed. 

Before running the pipeline, one may want to test it with a small vcf file (e.g. first 100,000 lines with CDS regions) to make sure things go alright. 




---------------------------------------------------------------------------------------------------------------------------------------------------------------------
An example of running

perl main.pl config.txt 
##### 1 configure file ############################
vcf_dir	/data/cjevol/test/input/
genome_dir	/data/cjevol/test/genome/
out_dir	/data/cjevol/test/output/
unfolded	1
Pop_tags	BWL
Outgr1	LP197014
Outgr2	LP197068
Outgr3	TMS17
num_inds	10
ploid	2
threads	4
filter	T
allSite	F
size	1000
write_genes	output_piNpiS_details
found reference genome /data/cjevol/test/genome/asm.cleaned.fasta.review.assembly.FINAL.fasta
found gff file /data/cjevol/test/genome/chr11.gff3

##### 2 convert vcf file ############################
perl util/convert-vcf2snp.pl /data/cjevol/test/genome/chr11.gff3 /data/cjevol/test/input/ tmp_snp 
Reading GFF
Start to convet VCF

Processing /data/cjevol/test/input/Quercus_all284_chr11_fltSNP.vcf.gz ...
Use of uninitialized value $fh in 1's complement (~) at util/convert-vcf2snp.pl line 47, <$fh> line 2172689. #### <--------- possible warning messages depending on system (no worries)
Done!

##### 3 Getting polymorphic sites ############################
perl calc-piNpiS-ithreads.pl -r /data/cjevol/test/genome/asm.cleaned.fasta.review.assembly.FINAL.fasta -g /data/cjevol/test/genome/chr11.gff3 -p snp -s 10 -d 2 -t 4 -o 1 -w focus_tmp 1> piNpiS_Ind10_diploid.txt 2> log_piNpiS_Ind10_diploid.txt

##### 4 Polarize frequency ############################
Get outgr aln
perl util/get_aln_from_snp.pl -r /data/cjevol/test/genome/asm.cleaned.fasta.review.assembly.FINAL.fasta -g /data/cjevol/test/genome/chr11.gff3 -p outgr_snp -s 3 -d 2 -t 4 -o 1 -w outgr_tmp 1>/dev/null 2>&1 
# prepare input (take a while) 
LP197014	1                  <----- can have more than one individuals (only consensus will be used for estSFS input. same for others) e.g. LP197014,LP197015,LP197016
LP197068	2
TMS17	3
Reading focus alignments focus_tmp/align/    <--------- this step takes quite a while 
..
Reading outgr alignments outgr_tmp/align/ and printing
Done
# running estSFS 
util/est-sfs version 2.03
Run   1 ML -24568.255547 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0806 r6[1] 0.0712 r6[2] 0.3425 r6[3] 0.3931 r6[4] 0.0648 r6[5] 0.0478 
Run   2 ML -24568.293031 k0 0.0348 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0926 
  r6[0] 0.0802 r6[1] 0.0705 r6[2] 0.3416 r6[3] 0.3958 r6[4] 0.0642 r6[5] 0.0477 
Run   3 ML -24568.350376 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0155 k4 0.0927 
  r6[0] 0.0828 r6[1] 0.0715 r6[2] 0.3407 r6[3] 0.3931 r6[4] 0.0646 r6[5] 0.0472 
Run   4 ML -25186.817212 k0 0.1317 k1 0.0234 k2 0.0004 k3 0.0175 k4 0.1431 
  r6[0] 0.0411 r6[1] 0.0490 r6[2] 0.2932 r6[3] 0.5586 r6[4] 0.0512 r6[5] 0.0068 
Run   5 ML -24568.259430 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0805 r6[1] 0.0712 r6[2] 0.3421 r6[3] 0.3939 r6[4] 0.0646 r6[5] 0.0476 
Run   6 ML -24568.255383 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0806 r6[1] 0.0712 r6[2] 0.3427 r6[3] 0.3930 r6[4] 0.0648 r6[5] 0.0478 
Run   7 ML -24568.261994 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0807 r6[1] 0.0716 r6[2] 0.3423 r6[3] 0.3932 r6[4] 0.0646 r6[5] 0.0476 
Run   8 ML -24794.153691 k0 0.0306 k1 0.0207 k2 0.0000 k3 0.0140 k4 0.0996 
  r6[0] 0.2545 r6[1] 0.0424 r6[2] 0.2969 r6[3] 0.3190 r6[4] 0.0418 r6[5] 0.0553 
Run   9 ML -24568.426310 k0 0.0346 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0837 r6[1] 0.0719 r6[2] 0.3416 r6[3] 0.3910 r6[4] 0.0646 r6[5] 0.0472 
Run  10 ML -26502.972162 k0 0.1560 k1 0.0811 k2 0.0000 k3 0.0114 k4 0.2463 
  r6[0] 0.0529 r6[1] 0.0357 r6[2] 0.5133 r6[3] 0.3512 r6[4] 0.0445 r6[5] 0.0124 
Overall ML -26502.9722 k0 0.0347 k1 0.0203 k2 0.0000 k3 0.0156 k4 0.0927 
  r6[0] 0.0806 r6[1] 0.0712 r6[2] 0.3427 r6[3] 0.3930 r6[4] 0.0648 r6[5] 0.0478 
20007 SNPs in total and 503 failed 

##### 5 Calculate piN/piS ############################
Done!
#################################








