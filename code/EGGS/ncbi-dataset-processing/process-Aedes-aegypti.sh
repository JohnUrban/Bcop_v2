#!/bin/bash
source EGGS-bash-utility-functions.txt 



## VAR
GCF=GCF_002204515.2
OUT=yellow_fever_mosquito.fasta

## AUTO FILL
DIR=ncbi_dataset/data
GTF=${DIR}/${GCF}/genomic.gtf
PROT=${DIR}/${GCF}/protein.faa


## EXECUTE

ncbi_pipeline ${GTF} ${PROT} ${OUT}

get_gene_bed_OCL_input ${GTF} >  yellow_fever_mosquito.OCLinput.tsv

get_seq_sizes ${DIR}/${GCF}/*.fna

get_chr_to_acc

