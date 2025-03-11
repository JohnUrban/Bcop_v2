
hgt=fungi ;
BEDDIR=../../input/
BEDPRE=closest-conserved-fly-gene-of-similar-size-to-alien-p450.CYPomeNames.3.sorted.genes.
BEDSFX=-only.bed
INBED=${BEDDIR}${BEDPRE}${hgt}${BEDSFX}

i=0 ;
while read line ; do
  let i++ ;
  idx=0${i} ;
  ARR=( echo $line ) ;
  BED=$( echo ${idx: -2}-${ARR[4]}-${hgt}.bed ) ;
  echo $BED ;
  echo $line | awk '{gsub(" ","\t"); print}' > ${BED} ;
done < ${INBED}
