import os, sys
from collections import defaultdict
import numpy as np
import pandas as pd

version='1.0.0'

help_string='''
Current Version: %s

Author: John M. Urban.

version 1.0
- Needed functions for reproducing Bcop_v2 paper EGGS analyses were copied from John Urban's orthoChainLinkerUtilities.py version 0.4.0.


''' % (version)





def read_size_file(fpath):
    '''fpath to 2-column tab-sep file with chrnames and lengths. i.e. what BEDtools calls a genome file.'''
    d = {}
    with open(fpath) as fh:
        for line in fh:
            line = line.strip().split()
            d[line[0]] = int(line[1])
    return d


def get_gene2og(fpath="Orthogroups.tsv"):
    og2genes = {}
    gene2og = {}
    i=0
    with open(fpath) as ogfh:
        for line in ogfh:
    #         print(line)
            i+=1
            line = line.strip('\n').split('\t')
            values = []
            for i in range(1, len(line)):
                for e in line[i].split(', '):
                    if e:
                        values.append(e)
            og2genes[line[0]] = [e for e in values if e ]
            vset = set()
            for e in values:
                if e:
                    vset.add(e)
                    vset.add(e.split('-')[0])
            for e in list(vset):
                try:
                    gene2og[e].append( line[0] )
                except:
                    gene2og[e] = [line[0]]
    return gene2og, og2genes

        

def process_OCL_input(genes, transtab=None, gene2og=None, gene_table_names=['chr', 'start', 'end', 'gene', 'strand'], trans_table_names=["gene", "protein","length"]):
    ## READ IN GENE TABLE
    fh = pd.read_csv(genes, sep="\t", names=gene_table_names)
    
    if transtab is not None:
        ## READ IN TRANSLATION TABLE
        trans = pd.read_csv(transtab, sep="\t", names=trans_table_names)

        ## For name-trans table, only keep rows with genes that were represented in the OrthoGroups
        trans = trans[trans['protein'].isin(list(gene2og.keys()))].reset_index().drop('index',axis=1)

        ## For gene loc input BED, only keep rows that are still represented in the name-trans table.
        fh = fh[fh['gene'].isin(list(trans['gene']))].reset_index().drop('index',axis=1)

        ## Perform the name-translation, updating the gene loc BED
        fh['protein'] = fh['gene'].map(lambda x: ','.join(list(set(list(trans.loc[trans['gene']==x,'protein'])))) )
    
        ## Add orthogroup information based on protein name.
        fh['og'] = fh['protein'].map(lambda x: ','.join( list(set( [omg for ele in x.split(',') for omg in gene2og[ele] ] ))) )

        ## Count number of orthogroups per protein. Should be 1. Can be >1 when >1 isoform was used.
        fh['Ngroups'] = fh['protein'].map(lambda x: len( set( [omg for ele in x.split(',') for omg in gene2og[ele] ] ) ) )

    else:
        fh = fh[fh['gene'].isin(list(gene2og.keys()))]   
        fh['og'] = fh['gene'].map(lambda x: ','.join(gene2og[x]))
        fh['Ngroups'] = fh['gene'].map(lambda x: len(gene2og[x]))
        
        
    ## Only keep proteins with a single OG
    fh = fh.loc[fh['Ngroups'] == 1, :].drop('Ngroups',axis=1)

    ## Ensure sorted
    fh = fh.sort_values(by=['chr','start','end']).reset_index().drop('index',axis=1)
    
    ## Ensure column order
    columns = ['chr','start','end','gene','strand','og'] ###,'protein']
    for e in list(fh.columns):
        if e not in columns:
            columns.append(e)
    fh = fh[columns]

    ## Return dataframe: use chr_chain_maker() to get input to OCL from this dataframe.
    return fh



def get_simple_paf_format(A, B, A_g, B_g, A_id="A", B_id="B"):
    '''
    Inputs into this function:
    - A     = Species A updated OCL BED-like dataframe (chr, start, end, gene, strand, og, protein).
    - B     = Species B updated OCL BED-like dataframe (chr, start, end, gene, strand, og, protein).
    - A_g   = Species A chrom size dataframe (Input 4 df output from python pre-processing).
    - B_g   = Species B chrom size dataframe (Input 4 df output from python pre-processing).
    - A_id = Species A identifer (usually first letter of genus and 3 letters of species name -- e.g. bcop)
                - Defaults to "A"
    - B_id = Species B identifier (usually first letter of genus and 3 letters of species name -- e.g. dmel)
                - Defaults to "B"
    
    
    Assumes the following steps already taken for two species in an OrthoFinder analysis.
    - All BASH pre-processing to get:
        - Input 1: proteomes/OrthoFinder/Results_*date*/Orthogroups/Orthogroups_SingleCopyOrthologues.tsv from the .txt file of same name.
        - Input 2: "OCL input" BED-like format (chr, start, end, gene name, strand) for each species
        - Input 3: "gene_name prot_name prot_length file" -- for each species
        - Input 4: chrom size (genome) file for each species
    - All python pre-processing steps:
        - "Input 1" is read in to get_gene2og function to produce two dicts: gene2og and og2genes.
        - "Input 2 and 3" read in using process_OCL_input() function along with gene2og output from above line.
            - Gives A and B (inputs to this function).
        - "Input 4" read in using read_size_file() function.
            - For A_g and B_g (inputs to this function).

    - Note: 
        - Species A will be considered the "Target" or "Subject" species (columns 6-9)
        - Species B will be considered the "Query" species (columns 1-4)
        - The strand value in column 5 is arbitrarily set to "+" for all rows. 
            - Species-specific strand information is retained in columns after column 12.
        - Columns 10-12 are usually # matches, alignment block length, and MAPQ (or bitscore for BLAST).
            - Here they are just dummy values -- all set to "1" for every row.
            - They are numeric so as to avoid conflicts in downstream analyses (e.g. in LAVE).
        - The first 12 PAF-like.
        - The last 7 columns are not PAF-like as they lack the XX:x: prefix formatting.
    Final columns:
        - 1 = chr_species_B
        - 2 = chr_size__species_B
        - 3 = start_species_B
        - 4 = end_species_B
        - 5 = dummy_strand
        - 6 = chr_species_A
        - 7 = chr_size_species_A
        - 8 = start_species_A
        - 9 = end_species_A
        - 10 = dummy
        - 11 = dummy
        - 12 = dummy
        - 13 = orthogroup (og)
        - 14 = gene_species_A
        - 15 = gene_species_B
        - 16 = protein_species_A
        - 17 = protein_species_B
        - 18 = strand_species_A
        - 19 = strand_species_B
        '''
    A_og = A.sort_values(by='og').reset_index(drop=True).set_index(keys='og')
    B_og = B.sort_values(by='og').reset_index(drop=True).set_index(keys='og')
    AB = A_og.join(B_og, lsuffix="_"+A_id, rsuffix="_"+B_id).sort_values(by=['chr_'+A_id, 'start_'+A_id, 'end_'+A_id]).reset_index()
    AB['chr_len_'+A_id] = AB['chr_'+A_id].map(lambda x: A_g[x])
    AB['chr_len_'+B_id] = AB['chr_'+B_id].map(lambda x: B_g[x])
    AB['dummy_strand'] = AB['chr_'+B_id].map(lambda x: '+')
    AB['dummy_match'] = AB['chr_'+B_id].map(lambda x: 1)
    AB['dummy_alnlen'] = AB['chr_'+B_id].map(lambda x: 1)
    AB['dummy_mapq'] = AB['chr_'+B_id].map(lambda x: 1)
    paf = AB[['chr_'+B_id,'chr_len_'+B_id,'start_'+B_id,'end_'+B_id,'dummy_strand','chr_'+A_id,'chr_len_'+A_id,'start_'+A_id,'end_'+A_id,'dummy_match','dummy_alnlen','dummy_mapq','og','gene_'+A_id,'gene_'+B_id, 'protein_'+A_id,'protein_'+B_id,'strand_'+A_id,'strand_'+B_id]]
    return paf



def get_simple_paf_format_direct_from_input_files(input1, input2A, input2B, input3A, input3B, input4A, input4B, A_id="A", B_id="B"):
    '''
    
    
    Assumes the following steps already taken for two species in an OrthoFinder analysis.
    - All BASH pre-processing to get:
        - Input 1: proteomes/OrthoFinder/Results_*date*/Orthogroups/Orthogroups_SingleCopyOrthologues.tsv from the .txt file of same name.
        - Input 2: "OCL input" BED-like format (chr, start, end, gene name, strand) for each species
        - Input 3: "gene_name prot_name prot_length file" -- for each species
        - Input 4: chrom size (genome) file for each species


    - This wrapper does All python pre-processing steps:
        - I.   "Input 1" is read in to get_gene2og function to produce two dicts: gene2og and og2genes.
        - II.  "Input 2 and 3" read in using process_OCL_input() function along with gene2og output from above line.
                - Gives A and B (inputs to this function).
        - III. "Input 4" read in using read_size_file() function.
                - For A_g and B_g (inputs to this function).
            
    - It then runs get_simple_paf_format() with inputs created inside the wrapper here:
        - A     = Species A updated OCL BED-like dataframe (chr, start, end, gene, strand, og, protein).
        - B     = Species B updated OCL BED-like dataframe (chr, start, end, gene, strand, og, protein).
        - A_g   = Species A chrom size dataframe (Input 4 df output from python pre-processing).
        - B_g   = Species B chrom size dataframe (Input 4 df output from python pre-processing).
        - A_id  = Species A identifer (usually first letter of genus and 3 letters of species name -- e.g. bcop)
                    - Defaults to "A"
        - B_id  = Species B identifier (usually first letter of genus and 3 letters of species name -- e.g. dmel)
                    - Defaults to "B"

    - Notes on final PAF-like output: 
        - Species A will be considered the "Target" or "Subject" species (columns 6-9)
        - Species B will be considered the "Query" species (columns 1-4)
        - The strand value in column 5 is arbitrarily set to "+" for all rows. 
            - Species-specific strand information is retained in columns after column 12.
        - Columns 10-12 are usually # matches, alignment block length, and MAPQ (or bitscore for BLAST).
            - Here they are just dummy values -- all set to "1" for every row.
            - They are numeric so as to avoid conflicts in downstream analyses (e.g. in LAVE).
        - The first 12 PAF-like.
        - The last 7 columns are not PAF-like as they lack the XX:x: prefix formatting.
    
    - Final columns:
        - 1 = chr_species_B
        - 2 = chr_size__species_B
        - 3 = start_species_B
        - 4 = end_species_B
        - 5 = dummy_strand
        - 6 = chr_species_A
        - 7 = chr_size_species_A
        - 8 = start_species_A
        - 9 = end_species_A
        - 10 = dummy
        - 11 = dummy
        - 12 = dummy
        - 13 = orthogroup (og)
        - 14 = gene_species_A
        - 15 = gene_species_B
        - 16 = protein_species_A
        - 17 = protein_species_B
        - 18 = strand_species_A
        - 19 = strand_species_B
        '''

    ## Get Orthogroups and gene members; create dictionaries mapping OGs to genes and genes to OGs.
    gene2og, og2genes = get_gene2og(fpath=input1)

    A = process_OCL_input(genes = input2A,
                           transtab = input3A,
                           gene2og = gene2og)

    B = process_OCL_input(genes = input2B,
                           transtab = input3B,
                           gene2og = gene2og)
    A_g = read_size_file(input4A)
    B_g = read_size_file(input4B)
    return get_simple_paf_format(A=A, B=B, A_g=A_g, B_g=B_g, A_id=A_id, B_id=B_id)



def write_from_paf(sp1, sp2, paf, suffix=None):
    if suffix is None:
        name = sp1 + "_" + sp2 + ".paf"
    else:
        name = sp1 + "_" + sp2 + "-" + suffix + ".paf"
    paf.to_csv(name, sep="\t", header=False, index=False)


def write_simple_paf_format_direct_from_input_files(input1, input2A, input2B, input3A, input3B, input4A, input4B, A_id="A", B_id="B", suffix=None):
    '''See help for "get_simple_paf_format_direct_from_input_files()".
        - This simply also writes the PAF using A_id, B_id, and optional additional info in "suffix".
        - It does also return the PAF object.
    '''
    paf = get_simple_paf_format_direct_from_input_files(input1, input2A, input2B, input3A, input3B, input4A, input4B, A_id, B_id)
    write_from_paf(sp1, sp2, paf, suffix=None)
    return paf

  


















