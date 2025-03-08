

F=XpX-gDNA-control-rep1-phasemap-q0.extrashort.sorted.txt


function extractChr {
  F=${1}
  CHR=${2}
  SFX=$( echo $F | awk '{gsub(/\./,"\t"); print "."$(NF)}' )
  BASE=$( basename ${F} ${SFX} )
  date > ${CHR}.err
  echo "${CHR} started..." >> ${CHR}.err
  awk -v "CHR=${CHR}" '$1==CHR && $3==CHR' ${F} | gzip -c > ${BASE}.${CHR}${SFX}.gz
  echo "${CHR} finished..." >> ${CHR}.err
  date >> ${CHR}.err
}


for CHR in X II III IV ; do 
  date
  echo ${CHR}
  extractChr ${F} ${CHR} &
done
wait
