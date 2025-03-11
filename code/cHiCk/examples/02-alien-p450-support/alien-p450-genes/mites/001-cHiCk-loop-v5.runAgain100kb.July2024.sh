## Absolute PATH to cHiCk Script
RUN=~/Bcop_v2/code/cHiCk/cHiCk/pipelines/whole-gene-pipeline.sh
## Relative PATH from this dir: ../../../../cHiCk/pipelines/whole-gene-pipeline.sh
## Relative PATH from Bcop_v2 dir: Bcop_v2/code/cHiCk/cHiCk/pipelines/whole-gene-pipeline.sh


## LINKS - Absolute Path to BEDPE file created by the Phase Map Pipeline.
LINKS=/Path/to/XO-phasemap-q0.sameChr.bedpe

## G size - Absolute path to Bcop_v2 "genome file". Can use one given in Bcop_v2 dir.
G=~/Bcop_v2/other/Bcop_v2.0.fasta.genome

## PSEUDOCOUNT FOR LOG10 ONLY
PSEUDO=1

## MAINDIR
MAINDIR=${PWD}

## INDIR
INDIR=${MAINDIR}/input



## LOOP
for BED in ${INDIR}/*.bed ; do
  ## OUT
  OUTPRE=$( basename ${BED} .bed )
  mkdir -p ${OUTPRE}
  echo "... moving into ${OUTPRE} dir." 1>&2
  cd ${OUTPRE}
  ## SIZE of each end.....
  for BINSIZE in 100000 ; do
    echo ${OUTPRE} ${BINSIZE} 1>&2
    ## EXECUTE
    mkdir -p ${BINSIZE}
    echo "... moving into ${BINSIZE} dir." 1>&2
    cd ${BINSIZE}
    bash ${RUN} ${LINKS} ${BED} ${G} ${BINSIZE} ${OUTPRE} ${PSEUDO} &
    echo "... moving back to sub-main dir." 1>&2
    cd ../
  done
  cd ${MAINDIR}
done
wait


exit
