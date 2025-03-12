# EGGS
Entropy of Gene Group Shuffling.

## Code used for Bcop_v2 paper can be found here.
- Updates to EGGS will be hosted at https://github.com/JohnUrban/EGGS
- Updates will include automating the process for others to do their own EGGS barplots, clustermaps, tables, and more.


### STEP1: Download genome datasets from NCBI
The NCBI directory structure for each species looks something like:
- ncbi_dataset/data/GCF*/
	- Inside of which are files like:
		- genomic.gtf 
		- protein.faa 
		- genomic.gff 
	- In some cases "GCA*" is needed, not "GCF*".



### STEP2: Processing NCBI datasets
- See: 
– e.g. Get longest isoforms of each protein for each proteome.
- Examples of how to do this with Bash functions can be found in the EGGS sub-directory named "ncbi-dataset-processing".
	- See "EGGS-bash-utility-functions.txt" and "process-*.sh" scripts.
	- The bash functions depend on some of the python scripts in the "utilities" sub-directory of "Bcop_v2/code/".
		- tableFilter.py ; extractFastxEntries.py ; fastaFormatter.py ; fxSize.py 
- Examples of how to further process those files for GeneSpace are also given.
	- Those commands will depend on first running the "process-*.sh" command(s).
	- You may also need some of the same python scripts from "Bcop_v2/code/utilities" and another:
		- setOps.py (although this line can be commented out; it was just a QC).	



### STEP3: Run OrthoFinder on longest isoforms from all 5 species.
orthofinder -t 16 -a 16 -f proteomes/ 1> ortho.out 2> ortho.err
-f <dir>        Directory with proteome FASTA files.
 -t <int>        Number of parallel sequence search threads [Default = 16]
 -a <int>        Number of parallel analysis threads




### STEP4: Extract SCOs across all five species into TSV file.
- Uses grep.py (or grep)
- After OrthoFinder done, go to Orthogroups subdirectory and run the extraction command.
	- OrthoFinder/Results_*/Orthogroups
		- Notice the following two commands are inside "( )"
		- ( head -n 1 Orthogroups.tsv ; 
		    grep.py -p Orthogroups_SingleCopyOrthologues.txt -f Orthogroups.tsv -c 1 -C 1 ) > Orthogroups_SingleCopyOrthologues.tsv
		- "grep.py" can be found in the EGGS utilities sub-directory.




### STEP5: Run EGGS Processing
- The Jupyter Notebook in the EGGS sub-directory shows how to use the EGGS python library (EGGS.py therein) to compute expected and observed entropy and MinMax Normalized Entropy, as well as how to generate the plots in the paper.





## OTHER NOTES:
- All Python Utilities are also part of (and originally from): https://github.com/JohnUrban/sciaratools2
	- Any further developments to those utilities will be done there.
	- The versions here are fossilized to preserve the analyses.
