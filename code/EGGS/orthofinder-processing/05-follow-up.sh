## To do after the PAF generation with the "04" jupyter notebook step.
## Uses provided - bhyg-newnames.txt

## Give chr names to Phyg chromosomes
translateTable.py -i bcop_bhyg.paf -d bhyg-newnames.txt -k 1 -v 2 -c 1 --force > bcop_bhyg-newnames.paf 

## Make an arbitrary PAF got Bcop vs itself
awk 'OFS="\t" {print $6,$7,$8,$9,$5,$6,$7,$8,$9,$10,$11,$12}' bcop_bhyg-newnames.paf  > bcop_bcop.paf
