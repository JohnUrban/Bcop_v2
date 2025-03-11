source ../cHiCk-functions.deprecated.txt

#grablinksfrombedends contig_103-contains-II9A.bed 1000 ../../XO-phasemap-q0.sameChr.bedpe 1> contig_103-ends.bedpe 2> contig_103-ends.bedpe.err



BINSIZE=10000
grablinksfrombedends contig_103-contains-II9A.bed ${BINSIZE} ../../XO-phasemap-q0.sameChr.bedpe 2> contig_103-sep-ends.${BINSIZE}bp.bedpe.err


#F=contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.bedpe ; FBASE=$( basename ${F} .bedpe ) ; for WINDOWS in ../windows-*bed ; do BASE=$( basename $WINDOWS .bed ) ; echo $BASE ; awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' ${F} | coverageBed -counts -a ${WINDOWS} -b - > ${FBASE}.counts-in-${BASE}.bedGraph ; done

## with removing contigs from same contig
#F=contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.bedpe ; FBASE=$( basename ${F} .bedpe ) ; for WINDOWS in ../windows-*bed ; do BASE=$( basename $WINDOWS .bed ) ; echo $BASE ; awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' ${F} | intersectBed -v -a - -b contig_103-contains-II9A.bed | coverageBed -counts -a ${WINDOWS} -b - > ${FBASE}.counts-in-${BASE}.bedGraph ; done


## with removing only intervals over the "end" source
##left
F=contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.bedpe ; F2=contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.LEFT.bed ; FBASE=$( basename ${F} .bedpe ) ; for WINDOWS in ../windows-*bed ; do BASE=$( basename $WINDOWS .bed ) ; echo $BASE ; awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' ${F} | intersectBed -v -a - -b ${F2} | coverageBed -counts -a ${WINDOWS} -b - > ${FBASE}.counts-in-${BASE}.bedGraph ; done
##right
F=contig_103-contains-II9A_XO-phasemap-q0.sameChr.right-10000bp.bedpe ; F2=contig_103-contains-II9A_XO-phasemap-q0.sameChr.right-10000bp.RIGHT.bed ; FBASE=$( basename ${F} .bedpe ) ; for WINDOWS in ../windows-*bed ; do BASE=$( basename $WINDOWS .bed ) ; echo $BASE ; awk 'OFS="\t" {print $1,$2,$3"\n"$4,$5,$6}' ${F} | intersectBed -v -a - -b ${F2} | coverageBed -counts -a ${WINDOWS} -b - > ${FBASE}.counts-in-${BASE}.bedGraph ; donewindows-100kb



## NOTE / WARN :: For future recon -- I updated "grablinksfrombedends" in functions.txt to do the window counts after working it out here.. so "grablinksfrombedends" won't run correctly as written above now.
