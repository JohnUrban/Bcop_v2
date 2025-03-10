#!/bin/bash
source EGGS-bash-utility-functions.txt 

## BEFORE RUNNING:
## - Some of the processed files produced by the "01-process*.sh" scripts or NCBI dataset directories.
##      - Those files are named in below commands.
##      - That is not necessary, but one needs to point those commands to those files.
## - Once you set up this script to point at the files, you can comment out the "exit" command below (next line).


## NCBI Datasets -- set this variable to wherever the NCBI datasets are.
NCBI_DATASETS=

## GENE PROT PROTLEN TABLES -- UPDATED FUNCTION
## - I used updated code to recreate some files. 



## BCOP - Paths from Bcop_v2 Repo::  Bcop_v2/genes/ncbi-Bcop_v1_to_v2-geneset-liftOff.gff AND Bcop_v2/other/Bcop_v2.0.fasta.genome
get_gene2protein_with_lengths_transtab_from_ncbi_gtf_20231004 ../../../genes/ncbi-Bcop_v1_to_v2-geneset-liftOff.gff > bcop-gene-prot-protlen.tsv
get_gene_bed_OCL_input_20231004 ../../../genes/ncbi-Bcop_v1_to_v2-geneset-liftOff.gff > bcop-all.chrNames.bed
cp ../../../other/Bcop_v2.0.fasta.genome  bcop.genome


## DMEL
get_gene2protein_with_lengths_transtab_from_ncbi_gtf_20231004 ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/ncbi_dataset/data/GCF_000001215.4/genomic.gtf > dmel-gene-prot-protlen.tsv
get_gene_bed_OCL_input_20231004 ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/ncbi_dataset/data/GCF_000001215.4/genomic.gtf > dmel-all.bed
##grep.py -f ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/seq.sizes.tsv -p ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/chr2acc.tsv --column 1 --patterns_column 2 > dmel-acc.genome
translateTable.py -i dmel-all.bed -c 1 -d ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/chr2acc.tsv -k 2 -v 1 --force > dmel-all.chrNames.bed 
##translateTable.py -i dmel-acc.genome -c 1 -d ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/chr2acc.tsv -k 2 -v 1 --force > dmel.genome
translateTable.py -i ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/seq.sizes.tsv -c 1 -d ${NCBI_DATASETS}/Drosophila_melanogaster-Release6_ncbi/chr2acc.tsv -k 2 -v 1 --force > dmel.genome

## AEDES
get_gene2protein_with_lengths_transtab_from_ncbi_gtf_20231004 ${NCBI_DATASETS}/Aedes_aegypti/ncbi_dataset/data/GCF_002204515.2/genomic.gtf > aedes-gene-prot-protlen.tsv
get_gene_bed_OCL_input_20231004 ${NCBI_DATASETS}/Aedes_aegypti/ncbi_dataset/data/GCF_002204515.2/genomic.gtf > aedes-all.bed
##grep.py -f ${NCBI_DATASETS}/Aedes_aegypti/seq.sizes.tsv -p ${NCBI_DATASETS}/Aedes_aegypti/chr2acc.tsv --column 1 --patterns_column 2 > aedes-acc.genome
translateTable.py -i aedes-all.bed -c 1 -d ${NCBI_DATASETS}/Aedes_aegypti/chr2acc.tsv -k 2 -v 1 --force > aedes-all.chrNames.bed 
##translateTable.py -i aedes-acc.genome -c 1 -d ${NCBI_DATASETS}/Aedes_aegypti/chr2acc.tsv -k 2 -v 1 --force > aedes.genome
translateTable.py -i ${NCBI_DATASETS}/Aedes_aegypti/seq.sizes.tsv -c 1 -d ${NCBI_DATASETS}/Aedes_aegypti/chr2acc.tsv -k 2 -v 1 --force > aedes.genome

## ANOPH
get_gene2protein_with_lengths_transtab_from_ncbi_gtf_20231004 ${NCBI_DATASETS}/Anopheles_gambiae/ncbi_dataset/data/GCF_000005575.2/genomic.gtf > anoph-gene-prot-protlen.tsv
get_gene_bed_OCL_input_20231004 ${NCBI_DATASETS}/Anopheles_gambiae/ncbi_dataset/data/GCF_000005575.2/genomic.gtf > anoph-all.bed
##grep.py -f ${NCBI_DATASETS}/Anopheles_gambiae/seq.sizes.tsv -p ${NCBI_DATASETS}/Anopheles_gambiae/chr2acc.tsv --column 1 --patterns_column 2 > anoph-acc.genome
translateTable.py -i anoph-all.bed -c 1 -d ${NCBI_DATASETS}/Anopheles_gambiae/chr2acc.tsv -k 2 -v 1 --force > anoph-all.chrNames.bed 
##translateTable.py -i anoph-acc.genome -c 1 -d ${NCBI_DATASETS}/Anopheles_gambiae/chr2acc.tsv -k 2 -v 1 --force > anoph.genome
translateTable.py -i ${NCBI_DATASETS}/Anopheles_gambiae/seq.sizes.tsv -c 1 -d ${NCBI_DATASETS}/Anopheles_gambiae/chr2acc.tsv -k 2 -v 1 --force > anoph.genome

## BHYG
####get_gene2protein_with_lengths_transtab_from_ncbi_gtf_20231004 ${NCBI_DATASETS}/Pseudolycoriella_hygida/ncbi_dataset/data/GCA_029228625.1/genomic.gff > bhyg-gene-prot-protlen.tsv
cp ${NCBI_DATASETS}/Pseudolycoriella_hygida/gene2protein-with-protlengths-transtab.tsv bhyg-gene-prot-protlen.tsv 
get_gene_bed_OCL_input_20231004 ${NCBI_DATASETS}/Pseudolycoriella_hygida/ncbi_dataset/data/GCA_029228625.1/genomic.gff  > bhyg-all.bed
####translateTable.py -i bhyg-all.bed -c 1 -d ${NCBI_DATASETS}/Pseudolycoriella_hygida/chr2acc.tsv -k 2 -v 1 --force > bhyg-all.chrNames.bed 
translateTable.py -i bhyg-all.bed -c 1 -d ${NCBI_DATASETS}/Pseudolycoriella_hygida/chr2acc.tsv -k 2 -v 1 --force | awk '{sub(" ","_"); sub("chromosome_",""); print}' > bhyg-all.chrNames.bed 
###translateTable.py -i ${NCBI_DATASETS}/Pseudolycoriella_hygida/seq.sizes.tsv -c 1 -d ${NCBI_DATASETS}/Pseudolycoriella_hygida/chr2acc.tsv -k 2 -v 1 --force > bhyg.genome
translateTable.py -i ${NCBI_DATASETS}/Pseudolycoriella_hygida/seq.sizes.tsv -c 1 -d ${NCBI_DATASETS}/Pseudolycoriella_hygida/chr2acc.tsv -k 2 -v 1 --force | awk '{sub(" ","_"); sub("chromosome_",""); print}' > bhyg.genome
