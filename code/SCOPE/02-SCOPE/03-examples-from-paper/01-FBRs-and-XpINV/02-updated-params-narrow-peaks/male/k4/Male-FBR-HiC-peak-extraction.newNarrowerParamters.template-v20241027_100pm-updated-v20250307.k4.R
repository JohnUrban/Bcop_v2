##############################################################################################
### SCOPE - Scatter Cluster Of Paired Ends - Analysis Pipeline.
##############################################################################################
### For now cite the following pre-print (or follow-up publication associated with it):
###    Chromosome-scale scaffolding of the fungus gnat genome (Diptera: Bradysia coprophila)
###    John M. Urban, Susan A. Gerbi, Allan C. Spradling
###    bioRxiv 2022.11.03.515061; doi: https://doi.org/10.1101/2022.11.03.515061
##############################################################################################
### What you will need:
### - Paired-end reads (typically from Hi-C) mapped to the genome and processed with our phase-map pipeline, or something similar.
### - The end product of the mapping pipeline should include a text file (optionally gzipped) with:
###    - 1 paired mate map positions per line with the following minimum set of columns:
###    - chr1 pos1 chr2 pos2 mapq
###       - This is like Juicer's "extra short" format.
###    - The text file can have more columns, but when you read it into R, make sure it has those 5 properly named.
###    - SCOPE currently only does intra-chromosomal analysis, and assumes it is working with a single chromosome.
###       - Therefore, the text file should be further partitioned into multiple files:
###           - each corresponding to the self-self interactions of a given chromosome.
###           - If this is not done, then it should be done in R by the user, else the SCOPE results and plots will correspond to a mixture of all chromosoems/sequences.
###             - Thus, if needed to do in R, look into subsetting dataframes - and subset such that chr1==chr2.
### - You will also want to know the length of the chromosome (or contig) being analyzed to use as the tig_size parameter.
###     - Or to specify a "maxpair" < that length (maximum length between pos1 and pos2)
##############################################################################################

setwd("~/Google-Drives/dr.john.urban/Carnegie/Projects/Bcop_v2_release/hic/inversion-X-aln/laveAnalysis/Oct2024/02-updated-params-narrow-peaks/male/update-20250307/k4")
system("mkdir -p SCOPE-figures")
system("mkdir -p SCOPE-tables")
system("pwd")
system("ls")

## DEPENDENCIES
library(KernSmooth)
library(lattice)

## SCOPE FUNCTIONS
source("../../../Jurban-HiC-peak-extraction-functions-version_20250305.R")


##############################################################
## VARIABLES
##############################################################
Xlen <- 70507862
# IIlen <- 58343183
# IIIlen <- 71047972
# IVlen <- 97081274

tig_size <- Xlen


##############################################################
## READ IN HIC DATA
##############################################################
# xlinks.female <- read.extra.short.pe(fh="../../../../July2024/XP-phasemap-q0.extrashort.sorted.X-chromosome-only.txt.gz")
xlinks.male <- read.extra.short.pe(fh="../../../../July2024/XO-phasemap-q0.extrashort.sorted.X-chromosome-only.txt.gz")



##############################################################
## DUMMY VARIABLE FOR XLINKS TO MAKE THIS EASIER
##############################################################
# xlinks <- xlinks.female
xlinks <- xlinks.male

dim(xlinks)
head(xlinks)

############################################################################################################################
## PARAM PICKING 
############################################################################################################################

## By default, I use hierarchical clustering because it is deterministic.
## However, k-means usually gives the identical or highly similar results.
## Moreover, k-means seems more robust to differences between samples.
##   So, when keeping the quantile threshold for peak calling the same between samples, one sample might include more "noise" from weaker dots (when using "h").
##   "h" clustering does not necessarily always behave s.t. the brightest dots become the first k clusters when cutting the tree.
##   "k" clustering does seem to choose the brightest dots as the "k" clusters; probably b/c it is more a function of minimizing variance with all data points.
##   Sometimes it is necessary to adjust (e.g. increase) the quantile threshold for both h and k to behave the same.
##   For smaller values of k, often you will want to use higher quantile thresholds as you're targeting fewer peaks.

############################################################################################################################
## PARAM PICKING- FBR-SPECIFIC NOTES
############################################################################################################################
## I originally did this analysis with lots of smoothing, using big 250kb windows.
## That analysis was done in July 2024, and was the described in the first version of these analyses in the bcopv2 paper and figures.
## I may or may not change those parameters to this....
## That analysis can be found here: ../../../July2024/breakpoint-and-FBR-HiC-peak-extraction-UPDATE.R
## It also was replicated and explored here: ../../breakpoint-and-FBR-HiC-peak-extraction-UPDATE-OCT24.R
## I updated those parameters to use 5Mb as the cut-off for mate pair selection rather than 1 Mb and 5 MB cluster filtering step.
## When using 5 Mb to filter mate pairs, one can also specify K=4.
## The result is the 4 clusters desired are returned without need for cluster filtering.
## The visual is slightly different of course, as it has a wider diagonal "gap" (1Mb vs 5Mb).

## The original parameters that I used in July 2024 
#origParams <- list(MAPQ=10, bx=250000, k=12, gridsize=1024, minpair=1000e3, quantile_thresh=0.998, minRegionSeparation=5e6)

## The updated parameters where minpair is 5 Mb (instead of 1 Mb) and k is 4 (not 12)
## origParams.adj1 is the intended parameter set for this file!!!!!!!!!!!!!!!!!!!!!!!!!!
## MAIN Difference between "male" and "female" parameter sets of the same name is k=3 vs k=4; males do not have the X' breakpoints.

## By default, I use hierarchical clustering because it is deterministic.
## However, k-means usually gives the identical or highly similar results.
## Moreover, k-means seems more robust to differences between males and females.
##   Males do not have the breakpoint dots.
##   So, when keeping the quantile threshold for peak calling the same, males include more "noise" from weaker dots.
##   "h" clustering does not necessarily behave s.t. the brightest dots become the first k clusters when cutting the tree.
##   "k" clustering does seem to choose the brightesr dots as the "k" clusters; probably b/c it is more a function of minimizing variance with all data points.
##   To use "h" with males, one needs to set a slightly higher quantile threshold since we are expecting fewer peaks, and will want fewer grid cells isolated.
##   It was only necessary to increase the quantile for the original set of parameters from 0.998 (females) to 0.999 (males).
##      Although it seem 0.999 works for females too.
##   Similarly, 0.9999 works for the narrower set for both it seems.
##   The only reason to keep original at 0.998 for female for now then is to preserve what I did in the Bcopv2 paper. 
##      Future work can use 0.999.



############################################################################################################################
## JOHNS VARIOUS PARAMETER SETTINGS - CAN BE DELETED FOR CLEAN TEMPLATE, LEAVING JUST ONE GENERIC SETTING.
############################################################################################################################

# origParams.adj1.female <- list(MAPQ=10, bx=250000, k=4, gridsize=1024, minpair=5e6, maxpair=NA, quantile_thresh=0.999, minRegionSeparation=5e6, clust.method="h")
# origParams.adj1.male   <- list(MAPQ=10, bx=250000, k=3, gridsize=1024, minpair=5e6, maxpair=NA, quantile_thresh=0.999, minRegionSeparation=5e6, clust.method="h")
narrower.FBR.4.Params.female <- list(MAPQ=10, bx=50e3, k=4, gridsize=2048, minpair=5e6, maxpair=NA, quantile_thresh=0.9999, minRegionSeparation=5e6, clust.method="h") ## 2024-10-24-v1; same results as narrower.FBR.2.Params; but is less aggressive.
# narrower.FBR.4.Params.male   <- list(MAPQ=10, bx=50e3, k=3, gridsize=2048, minpair=5e6, maxpair=NA, quantile_thresh=0.9999, minRegionSeparation=5e6, clust.method="h") ## 2024-10-24-v1; same results as narrower.FBR.2.Params; but is less aggressive.

###############################################################################
## DUMMY VARIABLE FOR PARAMETER SETTINGS
###############################################################################
## This dummy variable is simply to make copy/pasting the code easier:

## PURPOSELY USING "female" PARAMS HERE SIMPLY TO CHANGE K to 4, WHICH HELPS IDENTIFY THE SAME FBR SUMMITS.
useParams <- narrower.FBR.4.Params.female	


useParams

###############################################################################
## PIPELINE
###############################################################################
system("mkdir -p SCOPE-figures")

## CREATE CONTACT.GRID.GAPPED (DIAGONAL GAP FROM MIN INTERACTION DISTANCE FILTERING)
contactGrid.gapped <- contactGrid.creation(hic = xlinks,
                                           MAPQ=useParams$MAPQ,
                                           tig_size = tig_size, 
                                           bx = useParams$bx, 
                                           gridsize = useParams$gridsize, 
                                           minpair = useParams$minpair,
                                           maxpair = useParams$maxpair)



## EXTRACT PEAK CELLS
peaks <- contactGrid.cell.peak.calling(contactGrid = contactGrid.gapped, 
                                       quantile_thresh = useParams$quantile_thresh)

## CLUSTER PEAK CELLS TO DEFINE "DOTS": method="h" for hierarchical. Use "k" for K-means. Results are usually the same.
clusters <- contactGrid.clusterPeakCells.CallRegions(peaks = peaks,
                                                     k = useParams$k,
                                                     method=useParams$clust.method) 
length(clusters$zSummit) ## Number of clusters BEFORE filtering


## UPDATE CLUSTERS TO ONLY CONSIST OF THOSE WITH MINIMUM REGION SEPARATION OR MORE.
## - THIS STEP IS NOT NECESSARY WITH THESE PARAMETERS, BUT HERE AS PART OF A TUTORIAL OF SORTS.
clusters <- contactGrid.ClusteredRegions.Filter(clusters = clusters, 
                                                minRegionSeparation = useParams$minRegionSeparation,
                                                updateClusters = TRUE)
length(clusters$zSummit) ## Number of clusters AFTER filtering


## SHOW INFORMATION ON INTERACTING LOCI DEFINED IN EACH CLUSTER
contactGrid.ClusteredRegions.Filter(clusters = clusters, 
                                    minRegionSeparation = useParams$minRegionSeparation)


## PLOT CONTACT.GRID.GAPPED AND SHOW CLUSTER SQUARES (CORNERS ON DOTS).
pdf("SCOPE-figures/01a-gapped-contactGrid-with-cluster-squares.minmaxnorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.gapped,
                                         clusters = clusters, 
                                         minRegionSeparation = useParams$minRegionSeparation, 
                                         transform.fxn = min.max.normalize, seglwd = 1, segcol = "green")
dev.off()




## NOW ALSO CREATE COMPLETE UNGAPPED CONTACT GRID
contactGrid.complete <- contactGrid.creation(hic = xlinks,
                                             MAPQ=useParams$MAPQ,
                                             tig_size = tig_size, 
                                             bx = useParams$bx, 
                                             gridsize = useParams$gridsize, 
                                             minpair = 0) #### this is 0 here to get all interactions

## EXTRACT INFORMATION FROM CLUSTER (DOT) REGIONS FROM COMPLETE UNGAPPED CONTACT GRID
clusters.complete <- return.ClusterInfo.From.Different.ContactGrid(peaks, contactGrid.complete, clusters)
## AGAIN, CLUSTER FILTERING NOT NEEDED HERE (RETURNS INPUT) BUT SHOWN ANYWAY.
clusters.complete <- contactGrid.ClusteredRegions.Filter(clusters = clusters.complete, 
                                                         minRegionSeparation = useParams$minRegionSeparation,
                                                         updateClusters = TRUE)


## PLOT CONTACT.GRID.COMPLETE AND SHOW CLUSTER SQUARES (CORNERS ON DOTS).
## Using LOG10 is necessary for the "Complete Grid" if you want to see long-range structure.
## PLOT LOG10 MEDIAN NORMALIZED VERSION CONTACT.GRID.COMPLETE AND SHOW CLUSTER SQUARES (CORNERS ON DOTS).
pdf("SCOPE-figures/02-complete-contactGrid-with-cluster-squares.log10mednorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.complete,
                                         clusters = clusters.complete, 
                                         minRegionSeparation = useParams$minRegionSeparation,, 
                                         transform.fxn = log10.median.normalize,
                                         coltrio=c("violet","darkblue","white","red","yellow"), 
                                         segcol = "green", seglwd = 0.5, sigdig = 2)
dev.off()









## UPDATE CONTACT.GRIDS (GAPPED AND COMPLETE) WITH VECTORIZED PAIRWISE INTERACTION DISTANCE AND FREQUENCY INFORMATION (IntDists and IntFreqs).
contactGrid.gapped   <- vectorizeIntFreqsAndIntDists(contactGrid.gapped)
contactGrid.complete <- vectorizeIntFreqsAndIntDists(contactGrid.complete)

## MAKE LABELS FOR PLOT LEGENDS BELOW: THIS WILL HELP ID WHICH COLORS GO TO WHICH CLUSTERS.
legend.labels <- makeLabelsForClusters(clusters.complete); legend.labels


## COMPARE INTERACTION FREQUENCIES OF CLUSTERS (DOTS) TO BACKGROUND INTERACTION FREQUENCIES :: LIMITED TO ONLY INTERACTIONS WITH THE SAME DISTANCE
## GAPPED GRID
df.gapped.identical.dist   <- compareClusterContactFreqSummitsToBackground(clusters, contactGrid.gapped, 
                                                                           distBufferFactor = contactGrid.gapped$grid.cell.length, 
                                                                           scaleByBufferFactor = FALSE, enforce.nonzero.intfreq = TRUE)
## COMPLETE GRID
df.complete.identical.dist <- compareClusterContactFreqSummitsToBackground(clusters.complete, contactGrid.complete, 
                                                                           distBufferFactor = contactGrid.gapped$grid.cell.length, 
                                                                           scaleByBufferFactor = FALSE)

## PLOT LOG10 FE OF RESULTS AND SHOW STATISTICS (GAPPED AND COMPLETE GRID SHOULD GIVE SAME OR VERY SIMILAR RESULTS)
## FROM COMPLETE UNGAPPED "COMPLETE GRID"
pdf("SCOPE-figures/03-complete-contactGrid.clusterSummitLog10FoldEnrichentOverFreqsOfIdenticalLengths.pdf", width = 7, height = 7)
FE.upper.limit <- 1e6 ## This is used to constrain the x-axis if there are outliers. Most of the signal is typically much less than 6 orders of magnitude enriched. 
## subset of df.complete.identical.dist$FE > 1e-300 is to get rid of -inf values.
df.complete.identical.dist.stats <- plotDensityOfInteractionRatios(df.complete.identical.dist[df.complete.identical.dist$FE > 1e-300 & df.complete.identical.dist$FE <= FE.upper.limit,], 
                                                                   transform.fxn = log10, 
                                                                   add.legend = TRUE, 
                                                                   legend.labels = legend.labels$label, 
                                                                   legend.loc="topright"); 
df.complete.identical.dist.stats 
dev.off()

## PLOT RAW FE AND SHOW STATISTICS (GAPPED AND COMPLETE GRID SHOULD GIVE SAME OR VERY SIMILAR RESULTS)
## FROM COMPLETE GRID
pdf("SCOPE-figures/04-complete-contactGrid.clusterSummitFoldEnrichentOverFreqsOfIdenticalLengths.FE_le_100.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.complete.identical.dist[df.complete.identical.dist$FE > 1e-300 & df.complete.identical.dist$FE<=100, ], 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, 
                               legend.labels = legend.labels$label, 
                               legend.loc="topright")
dev.off()


## COMPARE INTERACTION FREQUENCIES OF CLUSTERS (DOTS) TO BACKGROUND INTERACTION FREQUENCIES :: COMPARED TO ENTIRE GRID.
## ONLY WORTH DOING ON COMPLETE UNGAPPED GRID. THE BLANK DIAGONAL IN THE GAPPED GRID CAUSES SOME ARTIFACTS. RESULTS ARE OTHERWISE SIMILAR.
df.complete.entire.grid   <- compareClusterContactFreqSummitsToBackground(clusters.complete, contactGrid.complete, 
                                                                          distBufferFactor = 1e12, 
                                                                          scaleByBufferFactor = FALSE,
                                                                          enforce.nonzero.intfreq = TRUE)
## PLOT LOG10 OF RESULTS AND SHOW STATISTICS FROM COMPLETE GRID ONLY (OFTEN HIGHLY SIMILAR TO ABOVE PLOTS WHEN INTDIST LONG)
## FROM COMPLETE GRID
pdf("SCOPE-figures/05-complete-contactGrid.clusterSummitLog10FoldEnrichentOverAllFreqsInGrid.pdf", width = 7, height = 7)
df.complete.entire.grid.stats <- plotDensityOfInteractionRatios(df.complete.entire.grid[df.complete.entire.grid$FE < 1e6, ], 
                                                                transform.fxn = log10, 
                                                                add.legend = TRUE, 
                                                                legend.labels = legend.labels$label, 
                                                                legend.loc="topright"); df.complete.entire.grid.stats
dev.off()
## PLOT RAW FOLD ENRICHMENT 
pdf("SCOPE-figures/06-complete-contactGrid.clusterSummitFoldEnrichentOverAllFreqsInGrid.FE_le_100.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.complete.entire.grid[df.complete.entire.grid$FE<=100, ], 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, 
                               legend.labels = legend.labels$label, 
                               legend.loc="topright")
dev.off()


## COMPARE INTERACTION FREQUENCIES OF CLUSTERS (DOTS) TO SHORT RANGE INTERACTION FREQUENCIES :: SHORT RANGE / DOTS.
## ONLY WORTH DOING ON COMPLETE UNGAPPED GRID. THE BLANK DIAGONAL IN THE GAPPED GRID CAUSES SOME ARTIFACTS. RESULTS ARE OTHERWISE SIMILAR.


df.complete.short.range <- compareClusterContactFreqSummitsToShortRange(clusters.complete, contactGrid.complete, 
                                                                        maxGridCellDist=1, interpretAsNtDistance=FALSE, 
                                                                        enforce.nonzero.intfreq=FALSE)

## THIS SHOWS HOW MUCH HIGHER ENRICHED THE BIN WITH SHORTEST RANGE CONTACTS ARE OVER CLUSTER SUMMIT(S) OF INTEREST (SHORT/CLUST), NOT CLUST/SHORT.
pdf("SCOPE-figures/07-complete-contactGrid.shortRangeIntFreqsFoldEnrichmentOverClusterSummitIntFreqs.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.complete.short.range, 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, legend.labels = legend.labels$label, legend.loc="topright")
dev.off()



## SCATTER PLOT OF ALL INTERACTIONS::  DISTANCES VS FREQUENCIES.

## LOG10 FREQS AND DISTS (LOG LOG PLOT) - SCALED BY SHORT RANGE
## USE ONLY FREQS > 0; i.e. NONZERO.
## FILTERING TINY VALUES FURTHER IMPROVES PLOT BY COLLAPSING Y-AXIS A BIT.
## TAKE THE TOP 99.95% OF NON-ZERO VALUES (EXCLUDE BOTTOM 0.05% OF INT FREQS > 0)
pdf("SCOPE-figures/08a-complete-contactGrid.scatterPlot.IntDists-vs-IntFreqs.freqsScaledByMedShortRange.log10X.log10Y.downsample100k.pdf", width = 7, height = 7)
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = log10, 
                        transform.fxn.dists = log10, 
                        down.sample = 1e5, 
                        scaleFreqsByMedShortRange = TRUE, 
                        grid.col = rgb(0,0,0,0.05), clust.col = rgb(1,0,0,1), 
                        enforce.nonzero.intfreq = TRUE, 
                        enforce.min.intfreq.quantile=TRUE, 
                        min.intfreq.quantile=0.0005) ##grid.col = "black", clust.col = "red"
dev.off()

pdf("SCOPE-figures/08b-complete-contactGrid.scatterPlot.IntDists-vs-IntFreqs.freqsScaledByMedShortRange.log10X.log10Y.all-data.pdf", width = 7, height = 7)
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = log10, 
                        transform.fxn.dists = log10,
                        scaleFreqsByMedShortRange = TRUE, 
                        grid.col = rgb(0,0,0,0.01), clust.col = rgb(1,0,0,1), 
                        enforce.nonzero.intfreq = TRUE, 
                        enforce.min.intfreq.quantile=TRUE, 
                        min.intfreq.quantile=0.0005)
dev.off()




###############################################################################
################## WRITING TABLES
###############################################################################
system("mkdir -p SCOPE-tables")

## ENSURING USE OF RIGHT VARS.
clusters <- clusters 
contactGrid <- contactGrid.gapped

## CHR IS WHAT WILL BE WRITTEN IN THE FIRST POSITION (CHR) OF BED FILES AND OTHER SIMILAR FILES.
CHR <- "X"




####################################
## MOST USEFUL INFO IS SUMMITS
####################################
## NICE THING ABOUT THE SUMMIT POSITION IS THAT IT STAYS CONSTANT OR NEARLY CONSTANT ACROSS A WHOLE BUNCH OF PARAMETER SETTINGS.
##   IN CONTRAST, THE CLUSTER BOUNDARIES CAN DIFFER SIMPLY DUE TO CHOICE OF K. 
##   THAT IN TURN CHANGES THE MEDIAN POSITION OF THE CLUSTER.
##   IN SOME PARAMETER SETTINGS, THE MEDIAN POSITION IS THE SAME AS THE SUMMIT POSITION, BUT NOT GUARANTEED.

N <- length(clusters$hx)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BED FILE SPANNING ENTIRE REGION BETWEEN POS1 AND POS2
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
clustSummits <- data.frame(chr = rep(CHR,N), 
                           start = round(clusters$xSummit), 
                           end = round(clusters$ySummit + contactGrid$grid.cell.length),
                           summitScore = clusters$zSummit)
clustSummits

# Write out.
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interaction.summit.bed')
write.table(clustSummits, file = fh, sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)





## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BEDPE FILE THAT CAN BE USED IN IGV TO SHOW ARCS CONNECTING POS1 AND POS2 (OR JUICEBOX TO SHOW OFF-DIAGONAL SQUARES AROUND DOTS)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
clustSummit.contact.arcs.bedpe <- data.frame(chr1 = rep(CHR,N), start1 = round(clusters$xSummit), end1 = round(clusters$xSummit + contactGrid$grid.cell.length),
                                             chr2 = rep(CHR,N), start2 = round(clusters$ySummit), end2 = round(clusters$ySummit + contactGrid$grid.cell.length), 
                                             name=rep(".",N), score=rep(clusters$zSummit, N), 
                                             strand1=rep(".",N), strand2=rep(".",N), 
                                             color=rep("\"255,0,0\"",N))
clustSummit.contact.arcs.bedpe

## NO HEADER FOR IGV
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interactions.summit.arcs.IGV.bedpe')
write.table(clustSummit.contact.arcs.bedpe, file = fh, sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

## WITH HEADER FOR JUICEBOX
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interactions.summit.off-diagonal-dot-squares.Juicebox.bedpe')
write.table(clustSummit.contact.arcs.bedpe, file = fh, sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BEDPE FILE THAT CAN BE USED IN IGV TO SHOW ARCS DEFINING POS1 AND POS2 BED COORDS (OR JUICEBOX TO SHOW ON-DIAGONAL SQUARES AT POS1 and POS2)
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
clustSummit.contact.regions.bedpe <- rbind(data.frame(chr1 = rep(CHR,N), start1 = round(clusters$xSummit), end1 = round(clusters$xSummit + contactGrid$grid.cell.length),
                                                      chr2 = rep(CHR,N), start2 = round(clusters$xSummit), end2 = round(clusters$xSummit + contactGrid$grid.cell.length), 
                                                      name=rep(".",N), score=rep(clusters$zSummit, N), 
                                                      strand1=rep(".",N), strand2=rep(".",N), 
                                                      color=rep("\"255,0,0\"",N)),
                                           data.frame(chr1 = rep(CHR,N), start1 = round(clusters$ySummit), end1 = round(clusters$ySummit + contactGrid$grid.cell.length),
                                                      chr2 = rep(CHR,N), start2 = round(clusters$ySummit), end2 = round(clusters$ySummit + contactGrid$grid.cell.length), 
                                                      name=rep(".",N), score=rep(clusters$zSummit, N), 
                                                      strand1=rep(".",N), strand2=rep(".",N), 
                                                      color=rep("\"255,0,0\"",N)))
clustSummit.contact.regions.bedpe <- clustSummit.contact.regions.bedpe[order(clustSummit.contact.regions.bedpe$start1,clustSummit.contact.regions.bedpe$end1),]

## NO HEADER FOR IGV
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interactions.summit.regions.IGV.bedpe')
write.table(clustSummit.contact.regions.bedpe, file = fh, sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

## WITH HEADER FOR JUICEBOX
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interactions.summit.on-diagonal-pos-squares.JuiceBox.bedpe')
write.table(clustSummit.contact.regions.bedpe, file = fh, sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## BEDPE FILE THAT CAN BE USED IN JUICEBOX TO SHOW ARC SQUARES CONNECTING POS1 AND POS2 
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
clustSummit.contact.squares.bedpe <- data.frame(chr1 = rep(CHR,N), start1 = round(clusters$xSummit), end1 = round(clusters$ySummit + contactGrid$grid.cell.length),
                                                chr2 = rep(CHR,N), start2 = round(clusters$xSummit), end2 = round(clusters$ySummit + contactGrid$grid.cell.length), 
                                                name=rep(".",N), score=rep(clusters$zSummit, N), 
                                                strand1=rep(".",N), strand2=rep(".",N), 
                                                color=rep("\"255,0,0\"",N))
clustSummit.contact.squares.bedpe


## WITH HEADER FOR JUICEBOX
fh <- paste0('SCOPE-tables/chr',CHR,'-major-interactions.summit.arc-squares.Juicebox.bedpe')
write.table(clustSummit.contact.squares.bedpe, file = fh, sep="\t", col.names = TRUE, quote = FALSE, row.names = FALSE)



## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
## TABLE FILE DETAILING BED COORDS OF POS1 AND POS2, BUT NOT BEDPE FORMAT.
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## 
clustSummitsPE <- data.frame(chr = rep(CHR,N), 
                             start1 = round(clusters$xSummit),
                             end1 = round(clusters$xSummit + contactGrid$grid.cell.length),
                             start2 = round(clusters$ySummit),
                             end2 = round(clusters$ySummit + contactGrid$grid.cell.length),
                             summitScore = clusters$zSummit)
clustSummitsPE

fh <- paste0('SCOPE-tables/chr',CHR,'-major-interaction.summit.paired.tsv')
write.table(clustSummitsPE, file = fh, sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

## Can use Awk with above file like this: awk 'OFS="\t" {print $1,$2,$3,$1"-dot_"NR"\n"$1,$4,$5,$1"-dot_"NR}' if you want.


## SUMMIT SIZE CHECK
slen <- clustSummitsPE$end1 - clustSummitsPE$start1
slen
contactGrid$grid.cell.length
slen/contactGrid$grid.cell.length  ## if slen=grd.cell.length, the ratio still might be slightly offset from 1 due to rounding above for pieces that make up slen.





######################################################
## FEEL FREE TO COMPARE SUMMIT POSITIONS TO MEDIAN POSITIONS.
## INTERACTING LOCUS CENTER (MEDIAN) LINKS
##  NOTE ALSO -- FARTHER BELOW IN DOC, MESS AROUND WITH ENTIRE REGIONS DEFINED BY MIN AND MAX POSITIONS.
######################################################

## MEDIAN BED
clustMedianPeak <- data.frame(chr = rep(CHR,N), 
                              start = round(clusters$hx), 
                              end = round(clusters$hy),
                              length = round(clusters$hy) - round(clusters$hx))
clustMedianPeak

## MEDIAN PAIRED BED
clustMedianPeakPE <- data.frame(chr = rep(CHR,N), 
                                start1 = round(clusters$hx),
                                end1 = round(clusters$hx + contactGrid$grid.cell.length),
                                start2 = round(clusters$hy),
                                end2 = round(clusters$hy + contactGrid$grid.cell.length),
                                length1 = round(clusters$hx + contactGrid$grid.cell.length) - round(clusters$hx),
                                length2 = round(clusters$hy + contactGrid$grid.cell.length) - round(clusters$hy))
clustMedianPeakPE



#############################################################################################################################################
## END OF ANALYSIS PIPELINE
#############################################################################################################################################


##############################################################################################
## BELOW ARE ALTERNATIVE PLOTS AND TABLES TO EXPLORE BUT EITHER REDUNDANT OR NOT ALWAYS USEFUL.
##############################################################################################



########################
## OTHER EXAMPLES OF PLOTTING THE CONTACT GRIDS AND INTERACTION FREQUENCY DISTRIBUTIONS.
########################

## PLOT LOG10 MEDIAN NORMALIZED VERSION CONTACT.GRID.COMPLETE WITHOUT SQUARES
pdf("SCOPE-figures/02-complete-contactGrid-without-squares.log10mednorm.pdf", width = 7, height = 7)
plotContactGrid.Raw(contactGrid = contactGrid.complete, 
                    transform.fxn = log10.median.normalize, 
                    coltrio=c("violet","darkblue","white","red","yellow"), 
                    sigdig = 2)
dev.off()



## PLOTTING THE COMPLETE GRID WITH SQUARES DENOTING LOCATION OF LONG-RANGE PEAKS.
## IT LOOKS BAD B/C DIAGONAL SO HIGH.
pdf("SCOPE-figures/02a-complete-contactGrid-with-cluster-squares.minmaxnorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.complete,
                                         clusters = clusters.complete, 
                                         minRegionSeparation = useParams$minRegionSeparation,
                                         transform.fxn = min.max.normalize, seglwd = 0.5, segcol = "green")
dev.off()



## PLOT LOG10 OF RESULTS AND SHOW STATISTICS (SHOULD BE THE SAME OR HIGHLY SIMILAR AS "COMPLETE GRID")
## FROM GAPPED GRID -- AGAIN: THIS SHOULD BE (OR USUALLY IS) REDUNDANT WITH THE SAME PLOT FROM THE UNGAPPED "COMPLETE GRID".
pdf("SCOPE-figures/03a-gapped-contactGrid.clusterSummitLog10FoldEnrichentOverFreqsOfIdenticalLengths.pdf", width = 7, height = 7)
df.gapped.identical.dist.stats <- plotDensityOfInteractionRatios(df.gapped.identical.dist, 
                                                                 transform.fxn = log10, 
                                                                 add.legend = TRUE, 
                                                                 legend.labels = legend.labels$label, 
                                                                 legend.loc="topright"); 
df.gapped.identical.dist.stats 
dev.off()


## PLOT RAW FOLD ENRICHMENT AND SHOW STATISTICS (GAPPED AND COMPLETE GRID SHOULD GIVE SAME OR VERY SIMILAR RESULTS)
## FROM GAPPED GRID-- AGAIN: THIS SHOULD BE (OR USUALLY IS) REDUNDANT WITH THE SAME PLOT FROM THE UNGAPPED "COMPLETE GRID".
pdf("SCOPE-figures/04a-gapped-contactGrid.clusterSummitFoldEnrichentOverFreqsOfIdenticalLengths.FE_le_100.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.gapped.identical.dist[df.gapped.identical.dist$FE<=100,], 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, 
                               legend.labels = legend.labels$label, 
                               legend.loc="topright")
dev.off()




########################
## OTHER EXAMPLES OF PLOTTING INTERACTION DISTANCES VS INTERACTION FREQUENCIES
########################

## RAW PLOT
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = identity, 
                        transform.fxn.dists = identity, 
                        down.sample = 1e4, 
                        scaleFreqsByMedShortRange = FALSE, 
                        grid.col = "black", clust.col = "red")

## LOG10 FREQS ALONE
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = log10, 
                        transform.fxn.dists = identity, 
                        down.sample = 1e4, 
                        scaleFreqsByMedShortRange = FALSE, 
                        grid.col = "black", clust.col = "red")


## LOG10 DISTS ALONE
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = identity, 
                        transform.fxn.dists = log10, 
                        down.sample = 1e4, 
                        scaleFreqsByMedShortRange = FALSE, 
                        grid.col = "black", clust.col = "red")


## LOG10 FREQS AND DISTS (LOG LOG PLOT)
contactGrid.scatterPlot(contactGrid = contactGrid.complete, 
                        clusters = clusters.complete, 
                        tooShort = -1, 
                        transform.fxn.freqs = log10, 
                        transform.fxn.dists = log10, 
                        down.sample = 1e4, 
                        scaleFreqsByMedShortRange = FALSE, 
                        grid.col = "black", clust.col = "red")

######################




##############################################################################################################
### ALTERNATIVE TABLES
##############################################################################################################

## ENSURING USE OF RIGHT VARS.
clusters <- clusters 
contactGrid <- contactGrid.gapped

## CHR IS WHAT WILL BE WRITTEN IN THE FIRST POSITION (CHR) OF BED FILES AND OTHER SIMILAR FILES.
CHR <- "X"
CHR <- "II"
CHR <- "III"
CHR <- "IV"


## THIS TABLE IS NOT GUARANTEED TO BE USEFUL AS CLUSTER MIN/MAX BOUNDARIES ARE HIGHLY DEPENDENT ON PARAMETER CHOICES.
## ALL INTERACTING LOCUS INFORMATION AS BEDPE
N<-length(clusters$hxmin); N
hierBEDPE <- data.frame(chr1     = rep(CHR,N), 
                        start1   = round(clusters$hxmin), 
                        end1     = round(clusters$hxmax), 
                        chr2     = rep(CHR,N), 
                        start2   = round(clusters$hymin), 
                        end2     = round(clusters$hymax), 
                        name     = rep(".",N), 
                        score    = rep(".",N), 
                        strand   = rep(".",N), 
                        color    = rep("255,0,0",N),
                        chrmed   = rep(CHR,N), 
                        startmed = round(clusters$hx), 
                        endmed   = round(clusters$hy),
                        startlen = round(clusters$hxmax)-round(clusters$hxmin), 
                        endlen   = round(clusters$hymax)-round(clusters$hymin), 
                        svlen    = round(clusters$hy)-round(clusters$hx))
hierBEDPE








