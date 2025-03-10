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


## After processing, next is running OrthoFinder using the protein.faa files.

## It has been a while since I ran this, dear reader, so it is possible that actually these NCBI processing steps can be done independent of OrthoFinder, or after it.
