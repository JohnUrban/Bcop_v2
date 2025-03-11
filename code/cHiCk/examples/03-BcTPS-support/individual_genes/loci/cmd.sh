cat <( awk 'OFS="\t" {print $0,"\"204,0,0\""}' 2/final-BcTPS-models.BcTPS-2_XO-phasemap-q0.sameChr.bedpe ) <( awk 'OFS="\t" {print $0,"\"0,204,0\""}' 3/final-BcTPS-models.BcTPS-3_XO-phasemap-q0.sameChr.bedpe ) <( awk 'OFS="\t" {print $0,"\"0,0,204\""}' 5/final-BcTPS-models.BcTPS-5_XO-phasemap-q0.sameChr.bedpe ) | sort -k1,1 -k2,2n -k5,5n > final-BcTPS-models.BcTPS-2_3_5_XO-phasemap-q0.sameChr.bedpe


cat input/* | awk 'OFS="\t" { MID=int(($2+$3)/2) ; START=MID-125000 ; E=MID+125000 ; print $1":"START"-"E,$4}' > locus-coords-for-IGV.250kb.tsv
cat input/* | awk '$1=="II"' | sortBed -i - | mergeBed -i - -d 100000 | awk 'OFS="\t" { MID=int(($2+$3)/2) ; START=MID-125000 ; E=MID+125000 ; print $1":"START"-"E,$4}' > multigene-chromosomeII-locus-coords-for-IGV.250kb.tsv 
