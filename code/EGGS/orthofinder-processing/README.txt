
## RECALL:
	- Genome datasets were downloaded from NCBI.
	- Each species had its own directory, inside of which was a NCBI directory structure that looked something like:
	- ncbi_dataset/data/GCF*/
		- Inside of which are files like:
			- genomic.gtf 
			- protein.faa 
			- genomic.gff 
	- In some cases "GCA*" is needed, not "GCF*".


	## STEPS THAT WERE DONE INSIDE THE DIR "".
	### 01-process*.sh
	- Does some processing from NCBI datasets, some of which will be over-written by updated functions in next step.
	- In the future, I will trim this down to the bare minimum needed to replicate.
	- For now, I had to scavenge directories for code, and retrace steps.


	### 02-further*.sh
	- Reprocesses some of the NCBI dataset stuff.
	- Generates files needed in subsequent steps.


## MAKE SURE TO HAVE RUN ORTHOFINDER, AND THEN MADE THE FILE NAMED:
- proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.tsv
- This is described in EGGS README.md
	- Uses grep.py (or grep)
	- After OrthoFinder done, go to Orthogroups subdirectory and run the extraction command.
		- OrthoFinder/Results_*/Orthogroups
			- Notice the following two commands are inside "( )"
			- ( head -n 1 Orthogroups.tsv ; 
			    grep.py -p Orthogroups_SingleCopyOrthologues.txt -f Orthogroups.tsv -c 1 -C 1 ) > Orthogroups_SingleCopyOrthologues.tsv
			- "grep.py" can be found in the EGGS utilities sub-directory.



## NEW STEPS THAT DEPEND ON ORTHOFINDER
### 03-further*.sh
- Takes input from previous steps, generates SCO lists and BED files for each genome.


## 04-prepare-for-analysis-in-R.ipynb 
- Jupyter notebook to take in outputs of previous steps as input to make PAF files.
	- Each line gives the chromosomal addresses for one SCO in species A and species B.
	- This is input to dot plots in R, using Lave: https://github.com/JohnUrban/lave


## 05-follow-up.sh 
- Does two small follow-up steps
- Uses provided file, bhyg-newnames.txt
