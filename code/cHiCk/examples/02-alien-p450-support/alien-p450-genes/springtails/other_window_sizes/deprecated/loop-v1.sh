## Script
RUN=/mnt/sequence/jurban/software/deepseep/useVersion/deepseep/pipelines/whole-gene-pipeline.sh

## LINKS AND WINDOWS
LINKS=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/XO-phasemap-q0.sameChr.bedpe

## G size
G=/mnt/sequence/jurban/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome

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

  ## SIZE of each end.....
  for BINSIZE in 10000 1000000 ; do
    echo ${OUTPRE} ${BINSIZE} 1>&2
    ## EXECUTE
    mkdir -p ${BINSIZE}
    echo "... moving into ${BINSIZE} dir." 1>&2
    bash ${RUN} ${LINKS} ${BED} ${G} ${BINSIZE} ${OUTPRE} ${PSEUDO}
    echo "... moving back to main dir." 1>&2
    cd ${MAIN}
  done
done



exit


## GREP BED NAME
##BED=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/BcTPS-support/indiv/input/final-BcTPS-models.BcTPS-5.bed
