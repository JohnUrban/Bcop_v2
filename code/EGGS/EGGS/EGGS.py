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


########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
##############################################################################################################################################
''' UTILITIES AND DOT-PLOT-RELATED.'''
##############################################################################################################################################
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 



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

  



########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 
##############################################################################################################################################
''' ENTROPY-RELATED.'''
##############################################################################################################################################
########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## ########## 


class SpeciesFiles20231011(object):
    def __init__(self, sp):
        ''' sp    = a prefix, typically all lower-case; typically first letter of genus and first 3 letters from species; examples: bcop, bhyg, dmel.'''
        self.sp = sp
        self._get_up()   ## Makes the first letter of "sp" capital. That is it.
        self._get_bed()  ## Constructs bed file handle; retrieved by user with get_bed() or get_input2().
        self._get_txt()
        self._get_genome()
    def _get_up(self):
        self.up = self.sp[0].upper() + self.sp[1:]
    def _get_bed(self):
        self.bed = "../"+self.up+"-SCOs.bed"
    def _get_txt(self):
        self.txt = "../"+self.up+"-SCOs.txt"
    def _get_genome(self):
        self.genome = "../../shared_files/" + self.sp + ".genome"
    def get_bed(self):
        return self.bed
    def get_txt(self):
        return self.txt
    def get_genome(self):
        return self.genome


class SpeciesFiles(object):
    def __init__(self, sp, shared_dir):
        ''' sp         = String; Species prefix, typically all lower-case; typically first letter of genus and first 3 letters from species; examples: bcop, bhyg, dmel.
            shared_dir = String; A single shared directory where *-SCOs.bed, *-SCOs.txt and *.genome files can be found for all species being analyzed.
                         - E.g. See Bcop_v2/code/EGGS/orthofinder-processing/'''
        self.sp = sp
        self.shared_dir = shared_dir.rstrip("/")
        self._get_up()   ## Makes the first letter of "sp" capital. That is it.
        self._get_bed()  ## Constructs bed file handle; retrieved by user with get_bed() or get_input2().
        self._get_txt()
        self._get_genome()
        
    def _get_up(self):
        self.up = self.sp[0].upper() + self.sp[1:]
    def _get_bed(self):
        self.bed = self.shared_dir + "/" + self.up + "-SCOs.bed"
    def _get_txt(self):
        self.txt = self.shared_dir + "/"  + self.up +  "-SCOs.txt"
    def _get_genome(self):
        self.genome = self.shared_dir + "/" + self.sp + ".genome"
    def get_bed(self):
        return self.bed
    def get_txt(self):
        return self.txt
    def get_genome(self):
        return self.genome




class SimplePaf(object):
    def __init__(self, sp1, sp2, input1, sp1id="sp1", sp2id="sp2", sp1_chr_names=None, sp2_chr_names=None, pseudo=0,
                 strip_sp1_names = False, strip_sp2_names=False, pre=["chr"], sfx=["L","R"]):
        self.input1 = input1 ## "../proteomes/OrthoFinder/Results_Oct10/Orthogroups/Orthogroups_SingleCopyOrthologues.tsv"
        self.sp1 = sp1
        self.sp2 = sp2
        self.sp1id = sp1id
        self.sp2id = sp2id
        self.sp1_chr_names = sp1_chr_names
        self.sp2_chr_names = sp2_chr_names
        self.pseudo = pseudo
        self._get_inputs()
        self._get_paf()
        self._set_chr_names() ## needs to be below _get_paf()
##        if strip_sp1_names:
##            self.strip_sp1_chr_names(pre=pre, sfx=sfx, paf=True, chrlist=True)
##        if strip_sp2_names:
##            self.strip_sp2_chr_names(pre=pre, sfx=sfx, paf=True, chrlist=True)
##        self._make_table()
##        self._entropy_pipeline()

    ##############################################################################################################################
    ### PUBLIC
    ##############################################################################################################################
    def write_paf(self, sp1pre=None, sp2pre=None, suffix=None):
        ''' Note sp1pre defaults to using self.sp1id, which in turn defaults to just "sp1";
                so having the option here gives you another chance to give a more specific name.
            Same for s2pre.
            Ultimately instead of writing sp1_sp2.paf, it can be Xxxx_Yyyy.paf, etc.
            Alternatively or in addition, suffix allows you to provide another string to make it unique.'''
        if sp1pre is None:
            sp1pre = self.sp1id
        if sp2pre is None:
            sp2pre = self.sp2id
        write_from_paf(sp1=sp1pre, sp2=sp2pre, paf=self.orig_paf, suffix=suffix)
        
        

    def strip_sp1_chr_names(self, pre=["chr"], sfx=["L","R"], paf=True, chrlist=True):
        '''Strips stuff off front and back of names. Useful when names present as chrom arms for example.'''
        if paf:
            self.paf['chr_'+self.sp1id+'_orig'] = self.paf['chr_'+self.sp1id]
            self.paf['chr_'+self.sp1id] = self.paf['chr_'+self.sp1id].map(lambda x: strip_string(x, pre, sfx))
        if chrlist:
            self.sp1_chr_names = [strip_string(e, pre, sfx) for e in self.sp1_chr_names]

    def strip_sp2_chr_names(self, pre=["chr"],sfx=["L","R"], paf=True, chrlist=True):
        if paf:
            self.paf['chr_'+self.sp2id+'_orig'] = self.paf['chr_'+self.sp2id]
            self.paf['chr_'+self.sp2id] = self.paf['chr_'+self.sp2id].map(lambda x: strip_string(x, pre, sfx))
        if chrlist:
            self.sp2_chr_names = [strip_string(e, pre, sfx) for e in self.sp2_chr_names]

    def reset_chr_names(self, sp1=False, sp2=False):
        if sp1:
            self.sp1_chr_names = list(self.paf.chr_sp1.unique())
        if sp2:
            self.sp2_chr_names = list(self.paf.chr_sp2.unique())

    
    def update_chr_names(self, sp1=False, sp2=False, sp1_chr_names=None, sp2_chr_names=None, force_all_chr = False):
        if sp1 or sp1_chr_names is not None:
            if sp1_chr_names is None or force_all_chr:
                self.reset_chr_names(sp1=True)
            else:
                self.sp1_chr_names = sp1_chr_names
        if sp2 or sp2_chr_names is not None:
            if sp2_chr_names is None or force_all_chr:
                self.reset_chr_names(sp2=True)
            else:
                self.sp2_chr_names = sp2_chr_names



    def update_pseudo(self, pseudo=None):
        if pseudo is not None:
            self.pseudo = pseudo




    def get_pseudo(self):
        return self.pseudo

    def get_sp1_chr_names(self):
        return self.sp1_chr_names

    def get_sp2_chr_names(self):
        return self.sp2_chr_names


    def get_counts_table(self):
        return self.tab
    
    def get_joint_probabilities_table(self):
        return self.joint_prob_tab


    def get_joint_entropy(self):
        return self.joint_entropy

    def get_random_joint_entropy(self):
        return self.random_joint_entropy

    def get_joint_entropy_relative_to_max(self):
        return self.joint_entropy/self.random_joint_entropy

    def get_joint_entropy_minmax_normalized_joint_approach(self):
        '''This is a min-max normalization where the score = (X-min)/(max-min).

        Here, the "min entropy" is an estimate of the entropy from joint probabilities from nearly perfect conservation by taking the average entropy from self-self joint prob matrix computed with pseudo-counts.

        And the "max entropy" is the entropy when gene shuffling is completely random with respect to chromosomes in both species given their gene (SCO) counts on those chromosomes.
            It is computed from joint probabilities that are computed as the product of pairwise marginal probabilities ; this is done with matrix multiplication using the marginal prob vectors.
            Here, the probability that a SCO is one chr_i in sp1 and chr_j in sp2 is marg_p(chr_j)*marg_p(chr_i) -- the proportion of SCOs on chr_i in sp1 multiplied by the proportion of SCOs on chr_j in sp2.

        That estimation process only works if a pseudo-count is given to the analysis.

        TODO: either force-compute with default pseudo in future; or return None value for self.conserved_joint_entropy, and for this function.

        '''
        X   = self.joint_entropy
        MIN = self.conserved_joint_entropy
        MAX = self.random_joint_entropy
        return (X-MIN)/(MAX-MIN)
        #return (self.joint_entropy-self.conserved_joint_entropy)/(self.random_joint_entropy-self.conserved_joint_entropy)


    def get_joint_entropy_minmax_normalized_margin_approach(self):
        '''This is a min-max normalization where the score = (X-min)/(max-min).

        Here, the "min entropy" is the average entropy from species-specific marginal entropies computed on the marginal probabilities from the input counts table (so includes pseudo-counts).
            This should give the same or very similar result as the joint approach.

        And the "max entropy" is the entropy when gene shuffling is completely random with respect to chromosomes in both species given their gene (SCO) counts on those chromosomes.
            It is computed from joint probabilities that are computed as the product of pairwise marginal probabilities ; this is done with matrix multiplication using the marginal prob vectors.
            Here, the probability that a SCO is one chr_i in sp1 and chr_j in sp2 is marg_p(chr_j)*marg_p(chr_i) -- the proportion of SCOs on chr_i in sp1 multiplied by the proportion of SCOs on chr_j in sp2.
        

        '''
        X   = self.joint_entropy
        MIN = self.avg_marg_entropy
        MAX = self.random_joint_entropy
        return (X-MIN)/(MAX-MIN)


    def get_joint_entropy_minmax_normalized_margin_approach_perfect(self):
        '''This is a min-max normalization where the score = (X-min)/(max-min).

        Here, the "min entropy" is the average entropy from species-specific marginal entropies computed on the marginal probabilities from the input counts table after removing pseudo-counts.
            This can fail if there is a chromosome (contig) present with 0 SCOs for a given species.
            Assuming all chromosomes in chrlist have 1 or more SCOs, this gives the exact marginal probs and exact marginal entropy values.
            Assuming a small pseudocount was used, it should give a very similar results to the regular margin approach.
            A caveat is that despite being the "exact" entropy value for perfect conservation (i.e. computed w/o pseudocounts) it may not be the perfect "min" control in a context where pseudo-counts were used.
                - To reiterate though, it should be nearly the same value as the approach with pseudocounts, assuming they were sufficiently small.

        And the "max entropy" is the entropy when gene shuffling is completely random with respect to chromosomes in both species given their gene (SCO) counts on those chromosomes.
            It is computed from joint probabilities that are computed as the product of pairwise marginal probabilities ; this is done with matrix multiplication using the marginal prob vectors.
            Here, the probability that a SCO is one chr_i in sp1 and chr_j in sp2 is marg_p(chr_j)*marg_p(chr_i) -- the proportion of SCOs on chr_i in sp1 multiplied by the proportion of SCOs on chr_j in sp2.
        

        '''
        X   = self.joint_entropy
        MIN = self.avg_marg_entropy_perfect
        MAX = self.random_joint_entropy
        return (X-MIN)/(MAX-MIN)
    

    def get_joint_expected_surprisal_normalized_to_random(self):
        return (2**self.joint_entropy)/(2**self.random_joint_entropy)

    

    def reanalyze(self, sp1_chr_names=None, sp2_chr_names=None, pseudo=None, force_all_chr = False,
                  strip_sp1_names = False, strip_sp2_names=False, pre=["chr"], sfx=["L","R"]):

        ## Possible new chromosome lists
        self.update_chr_names(sp1 = sp1_chr_names is not None or force_all_chr,
                              sp2 = sp2_chr_names is not None or force_all_chr,
                              sp1_chr_names = sp1_chr_names,
                              sp2_chr_names = sp2_chr_names,
                              force_all_chr = force_all_chr)

        ## Possibly strip names
        if strip_sp1_names:
            self.strip_sp1_chr_names(pre=pre, sfx=sfx, paf=True, chrlist=True)
        if strip_sp2_names:
            self.strip_sp2_chr_names(pre=pre, sfx=sfx, paf=True, chrlist=True)


        ## Compute (possibly new) gate
        self.table_gate = (self.paf.chr_sp1.map(lambda x: x in self.sp1_chr_names)) & (self.paf.chr_sp2.map(lambda x: x in self.sp2_chr_names))

        ## Compute new table
        self.tab = self.paf[self.table_gate].chr_sp1.groupby([self.paf.chr_sp1, self.paf.chr_sp2]).size().unstack().fillna(0)

        ## Add (possibly new) pseudocount
        self.update_pseudo(pseudo)
        self.tab += self.pseudo

        ## Re-run Entropy Pipeline
        self._entropy_pipeline()

        ## Individual steps replaced by Entropy Pipeline for now.
        ## Re-compute joint probs table (i.e. proportions of SCOs on X-and-Y).        
        ##self._get_joint_probabilities_from_table()
        ## Re-compute joint prob entropy
        ##self._compute_entropy_from_joint_probabilities()

    ##############################################################################################################################
    ### FOR __INIT__() ONLY
    ##############################################################################################################################
    def _get_inputs(self):
        ''' sp1 and sp2 are "SpeciesFiles" class objects.***
            *** at the moment I've used SpeciesFiles20231011, but I've set it up to be plug and play for a more generalized class.
        '''
        self.d = {}
        self.d['input1'] = self.input1
        self.d['input2A'] = self.sp1.get_bed()
        self.d['input2B'] = self.sp2.get_bed()
        self.d['input3A'] = self.sp1.get_txt()
        self.d['input3B'] = self.sp2.get_txt()
        self.d['input4A'] = self.sp1.get_genome()
        self.d['input4B'] = self.sp2.get_genome()



    def _get_paf(self): ## "sp1" and "sp2" purposely default instead of Xxxx identifiers.
        self.paf = get_simple_paf_format_direct_from_input_files(input1 = self.d['input1'],
                                                                 input2A = self.d['input2A'],
                                                                 input2B = self.d['input2B'],
                                                                 input3A = self.d['input3A'],
                                                                 input3B = self.d['input3B'],
                                                                 input4A = self.d['input4A'],
                                                                 input4B = self.d['input4B'],
                                                                 A_id=self.sp1id,
                                                                 B_id=self.sp2id)
        self.orig_paf = self.paf.copy() ## To allow other PAF to be manipulated.

    def _set_chr_names(self):
        self.update_chr_names(sp1=True, sp2=True, sp1_chr_names=self.sp1_chr_names, sp2_chr_names=self.sp2_chr_names, force_all_chr = False)


    def _make_table(self):
        self.table_gate = (self.paf.chr_sp1.map(lambda x: x in self.sp1_chr_names)) & (self.paf.chr_sp2.map(lambda x: x in self.sp2_chr_names))
        self.tab = self.paf[self.table_gate].chr_sp1.groupby([self.paf.chr_sp1, self.paf.chr_sp2]).size().unstack().fillna(0)
        self.tab += self.pseudo



    ## JOINT ENTROPY
    def _get_joint_probabilities_from_table(self):
        self.joint_prob_tab = get_joint_probabilities_from_table(tab=self.tab) #self.tab/self.tab.sum().sum()

    def _get_joint_prob_surpisal_table(self):
        self.surpisal_tab = get_surprisal_table_from_joint_probs(self.joint_prob_tab)

    def _get_weighted_joint_prob_surprisal_table(self):
        self.weight_surprisal_tab = get_weighted_joint_prob_surprisal_table(self.joint_prob_tab,
                                                                            self.surpisal_tab)
    def _compute_entropy_from_weighted_surprisals(self):
        self.joint_entropy = get_joint_prob_entropy(self.weight_surprisal_tab)

    def _compute_entropy_from_joint_probabilities(self):
        self.joint_entropy = get_entropy_from_joint_probabilities(prob=self.joint_prob_tab)
        ##-1 * ( self.joint_prob_tab * np.log2(self.joint_prob_tab) ).sum().sum()

    ## COMPLETELY RANDOM JOINT ENTROPY
    def _get_completely_random_joint_probabilities(self):
        self.random_joint_prob_tab = get_random_joint_probabilities_from_table(tab=self.tab)

    def _get_completely_random_joint_prob_surpisal_table(self):
        self.random_surpisal_tab = get_surprisal_table_from_joint_probs(self.random_joint_prob_tab)

    def _get_completely_random_weighted_joint_prob_surprisal_table(self):
        self.random_weight_surprisal_tab = get_weighted_joint_prob_surprisal_table(self.random_joint_prob_tab,
                                                                            self.random_surpisal_tab)

    def _compute_completely_random_entropy_from_weighted_surprisals(self):
        self.random_joint_entropy = get_joint_prob_entropy(self.random_weight_surprisal_tab)



    ### MARGIN ENTROPY
    def _get_margin_counts(self):
        self.sp1_marg_counts = get_margin_counts(self.tab,axis=1)
        self.sp2_marg_counts = get_margin_counts(self.tab,axis=0)

    def _get_marginal_probabilities_from_margin_counts(self):
        self.sp1_marg_probs = get_marginal_probabilities_from_margin_counts(self.sp1_marg_counts)
        self.sp2_marg_probs = get_marginal_probabilities_from_margin_counts(self.sp2_marg_counts)

    def _get_surprisal_table_from_marginal_probs(self):
        self.sp1_marg_surprisals = get_surprisal_table_from_marginal_probs(self.sp1_marg_probs)
        self.sp2_marg_surprisals = get_surprisal_table_from_marginal_probs(self.sp2_marg_probs)

    def _get_weighted_marginal_prob_surprisal_table(self):
        self.sp1_marg_weighted_surprisals = get_weighted_marginal_prob_surprisal_table(self.sp1_marg_probs,
                                                                                       self.sp1_marg_surprisals)
        self.sp2_marg_weighted_surprisals = get_weighted_marginal_prob_surprisal_table(self.sp2_marg_probs,
                                                                                       self.sp2_marg_surprisals)
    

    def _compute_marginal_entropy_from_weighted_surprisals(self):
        self.sp1_marg_entropy = get_marginal_prob_entropy(self.sp1_marg_weighted_surprisals)
        self.sp2_marg_entropy = get_marginal_prob_entropy(self.sp2_marg_weighted_surprisals)

    def _compute_marginal_entropy_from_marginal_probabilities(self):
        self.sp1_marg_entropy = get_entropy_from_marginal_probabilities(self.sp1_marg_probs)
        self.sp2_marg_entropy = get_entropy_from_marginal_probabilities(self.sp2_marg_probs)

    def _compute_average_marginal_entropy(self):
        self.avg_marg_entropy = np.mean([self.sp1_marg_entropy,
                                         self.sp2_marg_entropy])

    def _compute_perfect_conservation_entropy_scores_from_margins_without_pseudocounts(self):
        perfect_tab                   = self.tab - self.pseudo
        self.sp1_marg_entropy_perfect = get_marginal_prob_entropy_from_counts_table(tab = perfect_tab,
                                                                                    axis=1)
        self.sp2_marg_entropy_perfect = get_marginal_prob_entropy_from_counts_table(tab = perfect_tab,
                                                                                    axis=0)
        self.avg_marg_entropy_perfect = np.mean([self.sp1_marg_entropy_perfect,
                                                 self.sp2_marg_entropy_perfect])

    
    ## ESTIMATING AVG JOINT PROBABILITY ENTROPY FROM COMPLETELY CONSERVED SCENARIO. ((Note that entropy from joint probs converge on entropy from marg probs as conservation goes to 100% and off-diagonal cells go to 0 counts and 0 prob)).
    def _estimate_completely_conserved_entropy(self):
        self.self_self_dict = get_average_self_self_joint_entropy(tab=self.tab, pseudo=self.pseudo)
        self.conserved_joint_entropy = self.self_self_dict['H']

    
    ## JOINT ENTROPY PIPELINE
    def _entropy_pipeline(self):
        ## 2023-10-12
        ## In the absence of pseudo-counts, when ZEROS are encountered, it doesn't seem to give bad or wrong results here compared to using small pseudos....
        ## essentially the 0s that give log2 issues just cannot be computed and so do not have surprisals, for example, that added to the Entropy sum.
        ## And that is what the truth would converge on anyway -- as p in p*log2(p) goes to 0, the entire product approaches 0, thus nothing to add to the sum.
        ## While I see wisdom in addressing this issue anyway by anticipating it and having it do something on purpose (like add a tiny pseudo) or suppressing the "RuntimeWarning: divide by zero encountered in log2" warnings,
        ##  I don't see it as a priority for now.
        
        ## JOINT
        self._get_joint_probabilities_from_table()
        self._get_joint_prob_surpisal_table()
        self._get_weighted_joint_prob_surprisal_table()
        self._compute_entropy_from_weighted_surprisals()

        ## MAX ENTROPY GIVEN GENE DIST ACROSS CHROMS OF BOTH SPP
        self._get_completely_random_joint_probabilities()
        self._get_completely_random_joint_prob_surpisal_table()
        self._get_completely_random_weighted_joint_prob_surprisal_table()
        self._compute_completely_random_entropy_from_weighted_surprisals()

        ## MARG PROBS
        self._get_margin_counts()
        self._get_marginal_probabilities_from_margin_counts()
        self._get_surprisal_table_from_marginal_probs()
        self._get_weighted_marginal_prob_surprisal_table()
        self._compute_marginal_entropy_from_weighted_surprisals()
        self._compute_average_marginal_entropy()
        self._compute_perfect_conservation_entropy_scores_from_margins_without_pseudocounts()


        ## ESTIMATING LOWEST JOINT PROB ENTROPY (SHOULD BE SAME OR CONVERGE ON MARG PROB ENTROPY)
        self._estimate_completely_conserved_entropy()


    #def get_joint_probability_entropy_from_paf(paf, sp1names=None, sp2names=None, pseudo=0):
    #    tab = get_table_from_paf(paf, sp1names, sp2names,pseudo)
    #    prob = get_joint_probabilities_from_table(tab)
    #    return get_entropy_from_joint_probabilities(prob)
    





