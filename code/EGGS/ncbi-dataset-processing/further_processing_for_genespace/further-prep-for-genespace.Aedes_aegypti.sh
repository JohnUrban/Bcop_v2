BASE=yellow_fever_mosquito

TRANS=gene2protein-with-protlengths-transtab.longestOnly.tsv

TSV=${BASE}.protein-coding.OCLinput.tsv


# get lngest gene names
echo newest name strategy


grep ">" ${BASE}.fasta | awk '{sub(">",""); print $1}' > ${BASE}.fasta.names
grep.py -f gene2protein-with-protlengths-transtab.tsv -c 2 -p ${BASE}.fasta.names -C 1 > ${TRANS}
cut -f 1 ${TRANS} > ${BASE}.forGenespace.names.txt ## new strategy


# subset OCL input and translate names
## new strategy.
echo bed for genespace
grep.py -f ${BASE}.OCLinput.tsv -c 4 -p  ${BASE}.forGenespace.names.txt -C 1  | translateTable.py -i - -c 4 -d ${TRANS} -k 1 -v 2 | translateTable.py -i - -c 1 -d chr2acc.tsv -k 2 -v 1 --force > ${BASE}.forGenespace.bed



wc -l ${BASE}*
grep -c ">" ${BASE}*a


cut -f 4 ${BASE}.forGenespace.bed > del
grep ">" ${BASE}.fasta | awk '{sub(">",""); print $1}' > del2
wc -l del del2
setOps.py del del2 -i | wc -l
