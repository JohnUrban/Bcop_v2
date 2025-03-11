awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3, log($4+1)/log(10)}' contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.counts-in-windows-1kb.bedGraph > contig_103-contains-II9A_XO-phasemap-q0.sameChr.left-10000bp.counts-in-windows-1kb.log10.bedGraph
awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3, log($4+1)/log(10)}' contig_103-contains-II9A_XO-phasemap-q0.sameChr.right-10000bp.counts-in-windows-1kb.bedGraph > contig_103-contains-II9A_XO-phasemap-q0.sameChr.right-10000bp.counts-in-windows-1kb.log10.bedGraph


for BDG in *kb.bedGraph ; do BASE=$( basename ${BDG} .bedGraph ) ; echo $BASE ; awk 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3, log($4+1)/log(10)}' ${BDG} > ${BASE}.log10.bedGraph ; done
