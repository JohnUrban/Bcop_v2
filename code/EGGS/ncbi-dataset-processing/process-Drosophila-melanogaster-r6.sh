#!/bin/bash
source EGGS-bash-utility-functions.txt

## VAR
GCF=GCF_000001215.4
OUTPRE=fruit_fly_melanogaster

## AUTO FILL
OUT=${OUTPRE}.fasta
DIR=ncbi_dataset/data
GTF=${DIR}/${GCF}/genomic.gtf
PROT=${DIR}/${GCF}/protein.faa


## EXECUTE
#ncbi_pipeline ${GTF} ${PROT} ${OUT}
#get_gene_bed_OCL_input ${GTF} > ${OUTPRE}.OCLinput.tsv
#get_seq_sizes ${DIR}/${GCF}/*.fna
get_chr_to_acc



# get_gene_bed_OCL_input_20231004  ${GTF} > ${OUTPRE}.alt.OCLinput.tsv
#mv fruit_fly_melanogaster.OCLinput.tsv fruit_fly_melanogaster.old.OCLinput.tsv 
#mv fruit_fly_melanogaster.alt.OCLinput.tsv fruit_fly_melanogaster.OCLinput.tsv 


get_gene_bed_OCL_input_PROTEINCODING ${GTF} > ${OUTPRE}.protein-coding.OCLinput.tsv
