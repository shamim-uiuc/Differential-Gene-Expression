#!/usr/bin/bash

# Differentially expressed genes identification
# Output from tophat can be used to identify DEG genes

#cuffdiff –L tip,mid,base –T –b genome_file –c 50 –u  .gff file sample1_replicate1,sample1_replicate2 sample2_replicate1,sample2_replicate2

#-T Time series
#-L condition (here 3 condition tip,mid, base)
#-b bias correction
#-u Tells Cufflinks to do an initial estimation procedure to more accurately weight reads mapping to multiple locations in the genome
#-c The minimum number of alignments in a locus for needed to conduct significance testing
#library normalization method (default is geometric), There are other option quartile

Input_dir1=/shared_resources/shamim_share/RNAseq/tip_files
Input_dir2=/shared_resources/shamim_share/RNAseq/mid_files
Input_dir3=/shared_resources/shamim_share/RNAseq/base_files
gff_dir=/shared_resources/shamim_share/gff_files
genome_dir=/shared_resources/shamim_share/genomes/maize
Ouput_dir=/shared_resources/shamim_share/RNAseq/cuffdiff_output
	
cuffdiff -L tip,mid,base -T -c 50 -u -b ${genome_dir}/ZMGV3.fa -o ${Output_dir} ${gff_dir}/maize.gff ${Input_dir1}/accepted_hits.bam ${Input_dir2}/accepted_hits.bam ${Input_dir3}/accepted_hits.bam
