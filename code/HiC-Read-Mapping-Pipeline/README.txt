This mapping pipeline descends from recommendations from Phase Genomics on how to process Phase Genomics Hi-C data.
- Hence I named it "phase-map".


STEP 01:
This file is a clean template (but feel free to also see the examples below):
- 01phase-map.sh



STEP 02:
This is an example of how to separate the extrashort txt file by chromosome:
- partition-chromosomes.sh
The key is really a simple awk command:
- awk -v "CHR=${CHR}" '$1==CHR && $3==CHR' ${F}
This step can also be done in R, but doing it before loading into R means you can load smaller files.



EXAMPLES:
These files are the same except for what Hi-C data they pointed to.
- examples/phase-map.XO.sh
- examples/phase-map.XpX.sh
