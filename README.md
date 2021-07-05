# geneTree

For calculation of piN/piS ratio in geneTree project

1) plase notice that GTF or GFF file should be sorted by exon or CDS order (i.e. for negtiave strand the one with largest coordiante comes first)
2) most error comes from the gff format. If your CDS is not annotated by the format "Parent=([\w|\d|\.]+)" (which means only letters, numbers and dot) please modifiy line 36 in GenomeAnalysis.pm. 
contact me cjevol@zju.edu.cn
