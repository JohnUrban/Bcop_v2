#2023-10-10


## BEFORE RUNNING:
## - Some of the processed files produced by the "01-process*.sh" scripts were copied into a shared directory.
##	- Those files are named in below commands.
##	- That is not necessary, but one needs to point those commands to those files.
## - Once you set up this script to point at the files, you can comment out the "exit" command below (next line).
exit


SHARED=shared_files

####################################################################################
## PROTEIN TO GENE DICTIONARIES
####################################################################################


SCOs=proteomes/OrthoFinder/Results_Oct10/Orthogroups/Orthogroups_SingleCopyOrthologues.tsv
## Orthogroup	2=african_malaria_mosquito	3=black_fungus_gnat_Brazil	4=black_fungus_gnat_USA	5=fruit_fly_melanogaster	6=yellow_fever_mosquito

# Get Bcop SCOs
grep.py --file ${SHARED}/bcop-gene-prot-protlen.tsv --patterns ${SCOs} --column 2 --patterns_column 4 > Bcop-SCOs.txt


# Get Bhyg SCOs
grep.py --file ${SHARED}/bhyg-gene-prot-protlen.tsv --patterns ${SCOs} --column 2 --patterns_column 3 > Bhyg-SCOs.txt


# Get Dmel SCOs 
grep.py --file ${SHARED}/dmel-gene-prot-protlen.tsv --patterns ${SCOs} --column 2 --patterns_column 5 > Dmel-SCOs.txt


# Get Aedes SCOs
grep.py --file ${SHARED}/aedes-gene-prot-protlen.tsv --patterns ${SCOs} --column 2 --patterns_column 6 > Aedes-SCOs.txt


# Get Anoph SCOs
grep.py --file ${SHARED}/anoph-gene-prot-protlen.tsv --patterns ${SCOs} --column 2 --patterns_column 2 > Anoph-SCOs.txt


####################################################################################
## GENE LOC BED FILES
####################################################################################

## Bcop
grep.py --file ${SHARED}/bcop-all.chrNames.bed --patterns Bcop-SCOs.txt --column 4 --patterns_column 1 > Bcop-SCOs.bed

## Bhyg
grep.py --file ${SHARED}/bhyg-all.chrNames.bed --patterns Bhyg-SCOs.txt --column 4 --patterns_column 1 > Bhyg-SCOs.bed

# Dmel
grep.py --file ${SHARED}/dmel-all.chrNames.bed --patterns Dmel-SCOs.txt --column 4 --patterns_column 1 > Dmel-SCOs.bed

# Aedes
grep.py --file ${SHARED}/aedes-all.chrNames.bed --patterns Aedes-SCOs.txt --column 4 --patterns_column 1 > Aedes-SCOs.bed

# Anopheles
grep.py --file ${SHARED}/anoph-all.chrNames.bed --patterns Anoph-SCOs.txt --column 4 --patterns_column 1 > Anoph-SCOs.bed





##
translateTable.py -i bcop_bhyg.paf -d bhyg-newnames.txt -k 1 -v 2 -c 1 --force > bcop_bhyg-newnames.paf 
awk 'OFS="\t" {print $6,$7,$8,$9,$5,$6,$7,$8,$9,$10,$11,$12}' bcop_bhyg-newnames.paf  >bcop_bcop.paf
