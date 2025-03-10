#!/bin/bash
source EGGS-bash-utility-functions.txt


get_gene2poritein_with_lengths_from_ncbi_bhyg_GCA_gff  ncbi_dataset/data/GCA_029228625.1/genomic.gff 
select_longest_proteins_from_transtab
extract_longest_proteins_from_ncbi_prot_fa ncbi_dataset/data/GCA_029228625.1/protein.faa black_fungus_gnat_Brazil.fasta

get_gene_bed_OCL_input_20231004 ncbi_dataset/data/GCA_029228625.1/genomic.gff > bhyg.OCLinput.tsv
get_seq_sizes ncbi_dataset/data/GCA_029228625.1/*fna



get_chr_to_acc_bhyg ncbi_dataset/data/GCA_029228625.1/GCA_029228625.1_FCFRP_Bhyg_1.0_genomic.fna
