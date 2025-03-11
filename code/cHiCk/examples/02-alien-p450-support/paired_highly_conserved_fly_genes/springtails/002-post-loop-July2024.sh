## MAINDIR
MAINDIR=${PWD}

## INDIR
INDIR=${MAINDIR}/input



## LOOP
for BED in ${INDIR}/*.bed ; do
  ## OUT
  OUTPRE=$( basename ${BED} .bed )
  CHROM=$( head -n 1 ${BED} | awk '{print $1}' )
  mkdir -p "${CHROM}/allbw/100000/"
  echo "cp ${OUTPRE}/100000/*log10*bw ${CHROM}/allbw/100000/"
  cp ${OUTPRE}/100000/*log10*bw ${CHROM}/allbw/100000/
  #echo "mv ${OUTPRE} ${CHROM}"
done
