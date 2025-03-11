# Evidence for given genes
- Note some of the code in here might be older than the code found here: ../../cHiCk/ (( i.e. ~/Bcop_v2/code/cHiCk/cHiCk/ ))
- The code in this directory was used for the analysis in the paper, and eventually evolved into cHiCk.
- One can also try sourcing the newest cHiCk functions from there inside the bash scripts in this dir.

# Results
- Partial results (bigwigs) were left in the individual gene output directories.
- These are for convenience, but can be erased and re-generated with the code.


# Bash scripts
- input/cmd.sh 
	- Processes BED files to help set-up the input to cHiCk loop code.
	- Run this before 01-cHiCk-loop*.sh
- 001-cHiCk-loop-v5.runAgain100kb.July2024.sh
	- Runs the cHiCk code to get Hi-C interactions emanating from given gene(s).
- 002-post-loop-July2024.sh
	- Helps better organize some key files from 001-cHiCk-loop*.sh by chromosome.
