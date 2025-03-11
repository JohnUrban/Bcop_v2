
## hack for diff color links: https://github.com/igvteam/igv/issues/676
( echo -e "chr1\tx1\tx2\tchrom2\tstart2\tend2\tname\tscore\tstrand1\tstrand2\tcolor" ; cat <( awk 'OFS="\t" {print $0,"\"204,0,0\""}' 2/final-BcTPS-models.BcTPS-2_XO-phasemap-q0.sameChr.bedpe ) <( awk 'OFS="\t" {print $0,"\"0,204,0\""}' 3/final-BcTPS-models.BcTPS-3_XO-phasemap-q0.sameChr.bedpe ) <( awk 'OFS="\t" {print $0,"\"0,0,204\""}' 5/final-BcTPS-models.BcTPS-5_XO-phasemap-q0.sameChr.bedpe ) | sort -k1,1 -k2,2n -k5,5n ) > final-BcTPS-models.BcTPS-2_3_5_XO-phasemap-q0.sameChr.bedpe
