#!/bin/bash

## PARAMETERS MOST LIKELY IN NEED OF CHANGING
REF=	## FASTA OF GENOME

FAIDX=${REF}.fai  ## samtools faidx if needed

SIZES=${REF}.genome ## faSize -detailed ${REF} > ${REF}.sizes; if need be.

BWAIDX=${REF}  ## bwa index ref if needed --- KEEP AS REF IS SIDE-BY-SIDE IN SAME DIR; ELSE CHANGE!!!!

CPU=16

MINMAPQ=0 ## 0 allows downstream tools can handle mapq the way they want

OUTPRE=phasemap-output

CLEAN=true

## READS
R1=  ## /Path/to/Hic-reads.R1.fastq.gz
R2=  ## /Path/to/Hic-reads.R2.fastq.gz

## Software
BWA=		# ~/software/bwa/bwa/
SAMBLASTER=	# /mnt/sequence/jurban/software/samblaster/samblaster/
JAR=		# /mnt/sequence/jurban/software/yahs/juicer_tools.2.20.00.jar

## PATH
export PATH=${BWA}:${SAMBLASTER}:${PATH}


## OTHER FILES (NO NEED TO EDIT THIS)
BEDPE=${OUTPRE}.bedpe
XSHORT=${OUTPRE}.extrashort.sorted.txt
XSHORTHIC=${OUTPRE}.extrashort.sorted.hic


## PREP - Make FAI if need be
( if [ ! -f ${FAIDX} ]; then samtools faidx ${REF} ; fi ) > prep.err 2>&1

## PREP - Make sizes if need be
( if [ ! -f ${SIZES} ]; then faSize -detailed ${REF} > ${SIZES} ; fi ) >> prep.err 2>&1

## PREP - Make BWA idx if need be
( if [ ! -f ${BWAIDX}.bwt ]; then bwa index ${REF}  ; fi ) >> prep.err 2>&1

# Map, markdup, filter : 2136 removes UNMAP,MUNMAP,SECONDARY,SUPPLEMENTARY  (alts are 2060:UNMAP,MUNMAP,SUPPLEMENTARY ; and 268:UNMAP,MUNMAP,SECONDARY ; 2136 is catch-all)
(bwa mem -5SP -t 16 ${BWAIDX} ${R1} ${R2} | samblaster --removeDups | samtools view -S -h -b -F 2316 -q ${MINMAPQ} > aligned.bam ) 2> phase_map.err

# BEDtools BEDPE format
bedtools bamtobed -bedpe -i aligned.bam 1> ${BEDPE} 2>bedpe.err

## Juicer Extra Short Format
( awk 'BEGIN{OFS="\t"}{print $1, $2, $4, $5, $8}' ${BEDPE} | sort -k1,1 -k3,3 -k2,2n -k4,4n > ${XSHORT} ) > bedpe2xshort.err 2>&1

## Create .hic file
( java -jar -Xmx32G ${JAR} pre ${XSHORT} ${XSHORTHIC} ${SIZES} ) > juicerpre.err 2>&1

## Remove intermediate files
if $CLEAN ; then rm aligned.bam ${BEDPE} ${XSHORT} ; fi
