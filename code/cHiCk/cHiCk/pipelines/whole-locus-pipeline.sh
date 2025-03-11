## This is currently identical to the whole gene pipeline; just renamed to emphasize it does not need to be genes.

## HIC LINKS (BEDPE)
LINKS=${1}

## GREP BED NAME
BED=${2}

# bedtools-like genome file (chr size)
G=${3}

## SIZE of each end.....
BINSIZE=${4}


## OUT
OUTPRE=${5}

## PSEUDOCOUNT FOR LOG10 ONLY
PSEUDOCOUNT=1



## FUNCTIONS
SCRIPTDIR=$( dirname ${0} )
FXNS=${SCRIPTDIR}/../bash/cHiCk.functions.txt
BIN=${SCRIPTDIR}/../bin/
export PATH=${BIN}:${PATH}




source ${FXNS}

## v2 in use since circa Mar 13, 2024 ; descended from more raw code made circa Dec 2023/Jan 2024.
#grablinksfrombed_v2_cpp ${BED} ${BINSIZE} ${LINKS} ${G} ${PSEUDOCOUNT} 2> ${OUTPRE}.${BINSIZE}.bedpe.err

## Bumped to v3 July 23, 2024
grablinksfrombed_v3_cpp ${BED} ${BINSIZE} ${LINKS} ${G} ${PSEUDOCOUNT} 2> ${OUTPRE}.${BINSIZE}.bedpe.err

