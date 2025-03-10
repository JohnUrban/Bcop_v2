BASE=fruit_fly_melanogaster

TRANS=gene2protein-with-protlengths-transtab.longestOnly.tsv

#TSV=${BASE}.protein-coding.OCLinput.tsv
TSV=${BASE}.OCLinput.tsv


# get lngest gene names
echo newest name strategy


grep ">" ${BASE}.fasta | awk '{sub(">",""); print $1}' > ${BASE}.fasta.names
##grep.py -f gene2protein-with-protlengths-transtab.tsv -c 2 -p final_longest_protein_selections.txt -C 1 > ${TRANS}
grep.py -f gene2protein-with-protlengths-transtab.tsv -c 2 -p ${BASE}.fasta.names -C 1 > ${TRANS}
cut -f 1 ${TRANS} > ${BASE}.forGenespace.names.txt ## new strategy


# subset OCL input and translate names
## new strategy.
echo bed for genespace
grep.py -f ${TSV} -c 4 -p  ${BASE}.forGenespace.names.txt -C 1  | translateTable.py -i - -c 4 -d ${TRANS} -k 1 -v 2 | translateTable.py -i - -c 1 -d chr2acc.tsv -k 2 -v 1 --force > ${BASE}.forGenespace.bed



wc -l ${BASE}*
grep -c ">" ${BASE}*a


cut -f 4 ${BASE}.forGenespace.bed > del
grep ">" ${BASE}.fasta | awk '{sub(">",""); print $1}' > del2
wc -l del del2
setOps.py del del2 -i | wc -l


########

NPS="$( sort del | uniq -c | awk '$1>1 {print $2}' | paste -sd "|" - )"
grep -E ${NPS} fruit_fly_melanogaster.forGenespace.bed | sortBed -i - | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$3-$2}' | tableFilter.py -n 4 -s 6 | awk 'OFS="\t" {print $1,$2,$3,$4,$5}' > tmp.bed
grep -v -E ${NPS} fruit_fly_melanogaster.forGenespace.bed | cat - tmp.bed | sortBed -i - > tmp2.bed
mv tmp2.bed fruit_fly_melanogaster.forGenespace.bed
cut -f 4 fruit_fly_melanogaster.forGenespace.bed > fruit_fly_melanogaster.forGenespace.names.iter2.txt
extractFastxEntries.py -f ncbi_dataset/data/GCF_000001215.4/protein.faa -n fruit_fly_melanogaster.forGenespace.names.iter2.txt > fruit_fly_melanogaster.forGenespace.fasta


####
wc -l ${BASE}*
grep -c ">" ${BASE}*a


cut -f 4 ${BASE}.forGenespace.bed > del
grep ">" ${BASE}.fasta | awk '{sub(">",""); print $1}' > del2
wc -l del del2
setOps.py del del2 -i | wc -l



