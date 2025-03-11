## First run phase-map.sh pipeline to map Hi-C reads.
## Use the BEDPE output from that pipeline.
## Below that is "XO-phasemap-q0.bedpe".

exit

##############################################################################
## COMMANDS
## See results in: example-post-phase-map-results.tar.gz
##############################################################################


## Grab only links on the same chromosome (or contig) (intra-chromosomal)
awk '$1==$4' XO-phasemap-q0.bedpe > XO-phasemap-q0.sameChr.bedpe

## Grab only links from different chromosomes (inter-chromosomal)
awk '$1 ~ /^I|^X/ && $4 ~/X|II|III|IV/ && $1!=$4' XO-phasemap-q0.bedpe > XO-phasemap-q0.diffChr.bedpe


## Grab links connected to larger region around BcTPS on II
#II:28,132,369-28,288,970 -- broader
awk -v "A=28132369" -v "B=28288970" '$1==$4 && $1=="II" {if(  ( $2 >= A && $2 <= B ) || ( $3 >= A && $3 <= B ) || ( $5 >= A && $5 <= B ) || ( $6 >= A && $6 <= B) ){print }}' XO-phasemap-q0.sameChr.bedpe > BcTPS-on-II.broader.bedpe


## Grab links connected to larger region around BcTPS on II
#II:28,183,720-28,213,911 -- just around cluster of 3 genes
awk -v "A=28183720" -v "B=28213911" '$1==$4 && $1=="II" {if(  ( $2 >= A && $2 <= B ) || ( $3 >= A && $3 <= B ) || ( $5 >= A && $5 <= B ) || ( $6 >= A && $6 <= B) ){print }}' XO-phasemap-q0.sameChr.bedpe > BcTPS-on-II.bedpe


## Grab links connected to larger region around BcTPS on II
#III:69,318,749-69,334,122
awk -v "A=69318749" -v "B=69334122" '$1==$4 && $1=="III" {if(  ( $2 >= A && $2 <= B ) || ( $3 >= A && $3 <= B ) || ( $5 >= A && $5 <= B ) || ( $6 >= A && $6 <= B) ){print }}' XO-phasemap-q0.sameChr.bedpe > BcTPS-on-III.bedpe


## Grab links connected to larger region around BcTPS on II
#IV:87,227,363-87,235,222
awk -v "A=87227363" -v "B=87235222" '$1==$4 && $1=="IV" {if(  ( $2 >= A && $2 <= B ) || ( $3 >= A && $3 <= B ) || ( $5 >= A && $5 <= B ) || ( $6 >= A && $6 <= B) ){print }}' XO-phasemap-q0.sameChr.bedpe > BcTPS-on-IV.bedpe


## Combine all links from each chromosome above
#cat BcTPS-on-II.bedpe BcTPS-on-III.bedpe BcTPS-on-IV.bedpe > BcTPS-on-all-HiC-links.bedpe
#or
cat BcTPS-on-II.bedpe BcTPS-on-III.bedpe BcTPS-on-IV.bedpe | sortBed -i - > BcTPS-on-all-HiC-links.sorted.bedpe




