# geneTree

Please note if you have thousands of scaffolds and contigs in the genome, it is necessary to increase the number of openning file handles by running "ulimit -n 100000" in your linux or mac command line. 

update 2021-07:
1) one now can also specify genome file as '.fasta.gz', '.fa.gz' or '.fas.gz' and gff file with 'gff.gz', 'gtf.gz'.  












For calculation of piN/piS ratio in geneTree project

1) plase notice that GTF or GFF file should be sorted by exon or CDS order (i.e. for negtiave strand the one with largest coordiante comes first)
2) most error comes from the gff format. If your CDS is not annotated by the format "Parent=([\w|\d|\.]+)" (which means only letters, numbers and dot) please modifiy line 36 in GenomeAnalysis.pm. 
3) please refer to me orignial message for more details. 

Plase contact me cjevol@zju.edu.cn if help needed. 
