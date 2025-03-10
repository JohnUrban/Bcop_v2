#!/bin/bash
source EGGS-bash-utility-functions.txt

## VAR
GCF=GCF_000005575.2
OUTPRE=african_malaria_mosquito
OUT=african_malaria_mosquito.fasta

## AUTO FILL
DIR=ncbi_dataset/data
GTF=${DIR}/${GCF}/genomic.gtf
PROT=${DIR}/${GCF}/protein.faa


## EXECUTE
#ncbi_pipeline ${GTF} ${PROT} ${OUT}
#get_gene_bed_OCL_input ${GTF} > african_malaria_mosquito.OCLinput.tsv 
#get_seq_sizes ${DIR}/${GCF}/*.fna
#get_chr_to_acc


get_gene_bed_OCL_input_PROTEINCODING ${GTF} > ${OUTPRE}.protein-coding.OCLinput.tsv
