


RUN=/mnt/sequence/jurban/software/deepseep/useVersion/deepseep/pipelines/whole-gene-pipeline.sh

## LINKS AND WINDOWS
LINKS=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/XO-phasemap-q0.sameChr.bedpe

## GREP BED NAME
BED=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/BcTPS-support/indiv/input/final-BcTPS-models.BcTPS-4.bed

## G size
G=/mnt/sequence/jurban/base/data/sciara/phase/Bcop_v2.0/Bcop_v2.0_pkg/Bcop_v2.0.fasta.genome

## SIZE of each end.....
BINSIZE=10000

## OUT
OUTPRE=BcTPS-4


## PSEUDOCOUNT FOR LOG10 ONLY
PSEUDO=1

## EXECUTE

bash ${RUN} ${LINKS} ${BED} ${G} ${BINSIZE} ${OUTPRE} ${PSEUDO}




exit


