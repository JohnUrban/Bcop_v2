# EGGS
Entropy of Gene Group Shuffling.



# STEP1: Download genome datasets from NCBI
The NCBI directory structure for each species looks something like:
- ncbi_dataset/data/GCF*/
	- Inside of which are files like:
		- genomic.gtf 
		- protein.faa 
		- genomic.gff 
	- In some cases "GCA*" is needed, not "GCF*".


# STEP2: Processing NCBI datasets
– e.g. Get longest isoforms of each protein for each proteome.

	

# STEP3: Run OrthoFinder on longest isoforms from all 5 species.
orthofinder -t 16 -a 16 -f proteomes/ 1> ortho.out 2> ortho.err
-f <dir>        Directory with proteome FASTA files.
 -t <int>        Number of parallel sequence search threads [Default = 16]
 -a <int>        Number of parallel analysis threads

# STEP4: Extract SCOs across all five species into TSV file.
Uses grep.py (or grep)
- After OrthoFinder done, go to Orthogroups subdirectory and run the extraction command.
	- OrthoFinder/Results_*/Orthogroups
		- (head -n 1 Orthogroups.tsv ; grep.py -p Orthogroups_SingleCopyOrthologues.txt -f Orthogroups.tsv -c 1 -C 1 ) > Orthogroups_SingleCopyOrthologues.tsv

# STEP5: Run EGGS Processing
- In process of writing up documentation....


