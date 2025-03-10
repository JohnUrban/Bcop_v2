Genome datasets were downloaded from NCBI.

Each species had its own directory, inside of which was a NCBI directory structure that looked something like:
- ncbi_dataset/data/GCF*/
	- Inside of which are files like:
		- genomic.gtf 
		- protein.faa 
		- genomic.gff 
- In some cases "GCA*" is needed, not "GCF*".


## What each step/bash file does.

### 01-process*.sh
- Does some processing from NCBI datasets, some of which will be over-written by updated functions in next step.
- In the future, I will trim this down to the bare minimum needed to replicate.
	- For now, I had to scavenge directories for code, and retrace steps.


### 02-further*.sh
- Reprocesses some of the NCBI dataset stuff.
- Generates files needed in subsequent steps.

### 03-further*.sh
- Takes input from previous steps, generates SCO lists and BED files for each genome.


## 04-prepare-for-analysis-in-R.ipynb 
- Jupyter notebook to take in outputs of previous steps as input to make PAF files.
	- Each line gives the chromosomal addresses for one SCO in species A and species B.
	- This is input to dot plots in R, using Lave: https://github.com/JohnUrban/lave


## 05-follow-up.sh 
- Does two small follow-up steps
- Uses provided file, bhyg-newnames.txt
