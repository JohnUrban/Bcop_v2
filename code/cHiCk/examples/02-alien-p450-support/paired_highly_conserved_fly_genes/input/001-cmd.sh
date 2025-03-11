
# convert bed5 to bed6
awk 'OFS="\t" {print $1,$2,$3,$4,".",$5}' conserved-fly-genes.single-and-multi.bed > conserved-fly-genes.single-and-multi.bed6.bed 



# get closest conserved gene for each alien
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; closestBed -a ${BED}  -b conserved-fly-genes.single-and-multi.bed6.bed | awk 'OFS="\t" {print $7,$8,$9,$4"_"$10,$11,$12}' > closest-conserved-fly-gene-to-${BASE} ; done


# get closest gene that is of similar size within 1% diff
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk 'OFS="\t" {print $7,$8,$9,$4"_"$3-$2"."$10"_"$9-$8".dist_"$13,$11,$12}' ; done < ${BED} > closest-conserved-fly-gene-of-similar-size-to-${BASE} ; done 



## closest gene that is of similar size within 1% diff -- but print both the p450 and the conserved gene in sep files.... WITH COLORS
#black
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; CHROM=$( echo ${RECORD} | awk '{print $1}' ) ; mkdir -p ${CHROM} ; ( echo "track name="loci" description="FBRs and INV" visibility=2 itemRgb="On"" ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$2,$3,"0,0,0\n"$7,$8,$9,$10,$11,$12,$8,$9,"255,153,0"}' ) > ${CHROM}/${NAME}_${BEDLEN}.bed ; done < ${BED} ; done


## Above but using 50kb around center of gene instead of gene. 50 kb
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; CHROM=$( echo ${RECORD} | awk '{print $1}' ) ; mkdir -p ${CHROM} ; ( echo "track name="${NAME}" description="FBRs and INV" visibility=2 itemRgb="On"" ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk -v "HALFLEN=25000" 'OFS="\t" {M1=int((($2+$3)/2)+0.5) ; M2=int((($8+$9)/2)+0.5 ) ; print $1,M1-HALFLEN,M1+HALFLEN,$4,$5,$6,M1-HALFLEN,M1+HALFLEN,"0,0,0\n"$7,M2-HALFLEN,M2+HALFLEN,$10,$11,$12,M2-HALFLEN,M2+HALFLEN,"255,153,0"}' ) > ${CHROM}/${NAME}_${BEDLEN}.50kb.bed ; done < ${BED} ; done

# Above but using 100kb around center of gene instead of gene. 100 kb
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; CHROM=$( echo ${RECORD} | awk '{print $1}' ) ; mkdir -p ${CHROM} ; ( echo "track name="${NAME}" description="FBRs and INV" visibility=2 itemRgb="On"" ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk -v "HALFLEN=50000" 'OFS="\t" {M1=int((($2+$3)/2)+0.5) ; M2=int((($8+$9)/2)+0.5 ) ; print $1,M1-HALFLEN,M1+HALFLEN,$4,$5,$6,M1-HALFLEN,M1+HALFLEN,"0,0,0\n"$7,M2-HALFLEN,M2+HALFLEN,$10,$11,$12,M2-HALFLEN,M2+HALFLEN,"255,153,0"}' ) > ${CHROM}/${NAME}_${BEDLEN}.100kb.bed ; done < ${BED} ; done

# Above but using 150kb around center of gene instead of gene. 150 kb
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; CHROM=$( echo ${RECORD} | awk '{print $1}' ) ; mkdir -p ${CHROM} ; ( echo "track name="${NAME}" description="FBRs and INV" visibility=2 itemRgb="On"" ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk -v "HALFLEN=75000" 'OFS="\t" {M1=int((($2+$3)/2)+0.5) ; M2=int((($8+$9)/2)+0.5 ) ; print $1,M1-HALFLEN,M1+HALFLEN,$4,$5,$6,M1-HALFLEN,M1+HALFLEN,"0,0,0\n"$7,M2-HALFLEN,M2+HALFLEN,$10,$11,$12,M2-HALFLEN,M2+HALFLEN,"255,153,0"}' ) > ${CHROM}/${NAME}_${BEDLEN}.150kb.bed ; done < ${BED} ; done

# Above but using 200kb around center of gene instead of gene. 200 kb
for BED in ../../../p450-support/indiv/input/alien-p450.CYPomeNames.3.sorted.genes.*bed ; do BASE=$( basename ${BED} ) ; while read RECORD ; do NAME=$( echo ${RECORD} | awk '{print $4}' ) ; BEDLEN=$( echo ${RECORD} | awk '{print $3-$2}' ) ; CHROM=$( echo ${RECORD} | awk '{print $1}' ) ; mkdir -p ${CHROM} ; ( echo "track name="${NAME}" description="FBRs and INV" visibility=2 itemRgb="On"" ; echo -e "${RECORD}" | closestBed -d -a -  -b <( awk -v "BEDLEN=${BEDLEN}" '$3-$2 < 1.01*BEDLEN && $3-$2 > 0.99*BEDLEN' conserved-fly-genes.single-and-multi.bed6.bed ) | awk -v "HALFLEN=100000" 'OFS="\t" {M1=int((($2+$3)/2)+0.5) ; M2=int((($8+$9)/2)+0.5 ) ; print $1,M1-HALFLEN,M1+HALFLEN,$4,$5,$6,M1-HALFLEN,M1+HALFLEN,"0,0,0\n"$7,M2-HALFLEN,M2+HALFLEN,$10,$11,$12,M2-HALFLEN,M2+HALFLEN,"255,153,0"}' ) > ${CHROM}/${NAME}_${BEDLEN}.200kb.bed ; done < ${BED} ; done
