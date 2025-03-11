CONTIG=contig_103

## SIZE of each end.....
BINSIZE=10000

## GREP BED NAME
BEDBASE=${CONTIG}
BED=${BEDBASE}.bed


## SOURCE BED WHERE INTERVAL COMING FROM
SRCBED=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/bcopv1-contig-support/BcopV1_corrected_contigs_on_BcopV2.bed

## FUNCTIONS
FXNS=./cHiCk-functions.deprecated.txt

## LINKS AND WINDOWS
LINKS=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/XO-phasemap-q0.sameChr.bedpe
WINDOWS=/mnt/sequence/jurban/bcop_v2_paper/hic/analyses/Bcop_v2/XO/q0/rerun-20240122/bcopv1-contig-support/windows-1kb.bed




source ${FXNS}

grep -w ${CONTIG} ${SRCBED} > ${BED}

grablinksfrombed  ${BED} ${BINSIZE} ${LINKS} ${WINDOWS} 2> ${BEDBASE}-whole-interval.${BINSIZE}bp.bedpe.err
