hgt=springtails ; i=0 ; while read line ; do let i++ ; idx=0${i} ; ARR=( echo $line ) ; BED=$( echo ${idx: -2}-${ARR[4]}-${hgt}.bed ) ; echo $BED ; echo $line | awk '{gsub(" ","\t"); print}' > ${BED} ; done < ../../input/alien-p450.CYPomeNames.3.sorted.genes.springtails-only.bed 



