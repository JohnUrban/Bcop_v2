setwd("~/Google-Drives/dr.john.urban/Carnegie/Projects/Bcop_v2_release/hic/inversion-X-aln/laveAnalysis/Oct2024/01-orig-params-update/female/")

## TEMPLATE VERSION 2024-10-27-100pm.
## With modifications needed to get things to work, such as enforcing.nonzero, and limiting FE levels in some place..

## FOR DEVELOPMENT CODE THAT LEAD TO THE TIDY PIPELINE HEREIN, SEE ../

library(KernSmooth)
library(lattice)

## LOCAL FUNCTIONS
source("../../Jurban-HiC-peak-extraction-functions-version_20241023.R")



## ENV VARIABLES
Xlen <- 70507862
tig_size <- Xlen

## READ IN HIC DATA
xlinks.female <- read.extra.short.pe(fh="../../../July2024/XP-phasemap-q0.extrashort.sorted.X-chromosome-only.txt.gz")
# xlinks.male <- read.extra.short.pe(fh="../../../July2024/XO-phasemap-q0.extrashort.sorted.X-chromosome-only.txt.gz")

## DUMMY VARIABLE FOR XLINKS TO MAKE THIS EASIER
xlinks <- xlinks.female
# xlinks <- xlinks.male
dim(xlinks)


## PARAM PICKING
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
origParams.adj1.female <- list(MAPQ=10, bx=250000, k=4, gridsize=1024, minpair=5e6, quantile_thresh=0.999, minRegionSeparation=5e6, clust.method="h")
# origParams.adj1.male   <- list(MAPQ=10, bx=250000, k=3, gridsize=1024, minpair=5e6, quantile_thresh=0.999, minRegionSeparation=5e6, clust.method="h")
# narrower.FBR.4.Params.female <- list(MAPQ=10, bx=50e3, k=4, gridsize=2048, minpair=5e6, quantile_thresh=0.9999, minRegionSeparation=5e6, clust.method="h") ## 2024-10-24-v1; same results as narrower.FBR.2.Params; but is less aggressive.
# narrower.FBR.4.Params.male   <- list(MAPQ=10, bx=50e3, k=3, gridsize=2048, minpair=5e6, quantile_thresh=0.9999, minRegionSeparation=5e6, clust.method="h") ## 2024-10-24-v1; same results as narrower.FBR.2.Params; but is less aggressive.


## DUMMY VARIABLE FOR PARAMETERS
## This dummy variable is simply to make copy/pasting the code easier:
useParams <- origParams.adj1.female
# useParams <- narrower.FBR.4.Params.female
# useParams <- origParams.adj1.male
# useParams <- narrower.FBR.4.Params.male
useParams

###############################################################################
## PIPELINE
###############################################################################

## CREATE CONTACT.GRID.GAPPED (DIAGONAL GAP FROM MIN INTERACTION DISTANCE FILTERING)
contactGrid.gapped <- contactGrid.creation(hic = xlinks,
                                           MAPQ=useParams$MAPQ,
                                           tig_size = tig_size, 
                                           bx = useParams$bx, 
                                           gridsize = useParams$gridsize, 
                                           minpair = useParams$minpair)



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
pdf("01a-gapped-contactGrid-with-cluster-squares.minmaxnorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.gapped,
                                         clusters = clusters, 
                                         minRegionSeparation = useParams$minRegionSeparation, 
                                         transform.fxn = min.max.normalize)
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
## IT LOOKS BAD B/C DIAGONAL SO HIGH.
pdf("02a-complete-contactGrid-with-cluster-squares.minmaxnorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.complete,
                                         clusters = clusters.complete, 
                                         minRegionSeparation = useParams$minRegionSeparation,
                                         transform.fxn = min.max.normalize)
dev.off()


## PLOT LOG10 MEDIAN NORMALIZED VERSION CONTACT.GRID.COMPLETE AND SHOW CLUSTER SQUARES (CORNERS ON DOTS).
pdf("02b-complete-contactGrid-with-cluster-squares.log10mednorm.pdf", width = 7, height = 7)
plotContactGrid.WithClusteredPeakRegions(contactGrid = contactGrid.complete,
                                         clusters = clusters.complete, 
                                         minRegionSeparation = useParams$minRegionSeparation,, 
                                         transform.fxn = log10.median.normalize,
                                         coltrio=c("violet","darkblue","white","red","yellow"), 
                                         segcol = "green", seglwd = 1, sigdig = 2)
dev.off()



## PLOT LOG10 MEDIAN NORMALIZED VERSION CONTACT.GRID.COMPLETE WITHOUT SQUARES
pdf("02c-complete-contactGrid.log10mednorm.pdf", width = 7, height = 7)
plotContactGrid.Raw(contactGrid = contactGrid.complete, 
                    transform.fxn = log10.median.normalize, 
                    coltrio=c("violet","darkblue","white","red","yellow"), 
                    sigdig = 2)
dev.off()





## UPDATE CONTACT.GRIDS (GAPPED AND COMPLETE) WITH VECTORIZED PAIRWISE INTERACTION DISTANCE AND FREQUENCY INFORMATION (IntDists and IntFreqs).
contactGrid.gapped   <- vectorizeIntFreqsAndIntDists(contactGrid.gapped)
contactGrid.complete <- vectorizeIntFreqsAndIntDists(contactGrid.complete)

## MAKE LABELS FOR PLOT LEGENDS BELOW: THIS WILL HELP ID WHICH COLORS GO TO WHICH CLUSTERS.
legend.labels <- makeLabelsForClusters(clusters.complete); legend.labels


## COMPARE INTERACTION FREQUENCIES OF CLUSTERS (DOTS) TO BACKGROUND INTERACTION FREQUENCIES :: LIMITED TO ONLY INTERACTIONS WITH THE SAME DISTANCE
## GAPPED GRID
df.gapped.identical.dist   <- compareClusterContactFreqSummitsToBackground(clusters, contactGrid.gapped, 
                                                                           distBufferFactor = 0, 
                                                                           scaleByBufferFactor = FALSE, enforce.nonzero.intfreq = TRUE)
## COMPLETE GRID
df.complete.identical.dist <- compareClusterContactFreqSummitsToBackground(clusters.complete, contactGrid.complete, 
                                                                           distBufferFactor = 0, 
                                                                           scaleByBufferFactor = FALSE)

## PLOT LOG10 OF RESULTS AND SHOW STATISTICS FROM BOTH (SHOULD BE THE SAME OR HIGHLY SIMILAR)
## FROM GAPPED GRID
pdf("03a-gapped-contactGrid.clusterSummitLog10FoldEnrichentOverFreqsOfIdenticalLengths.pdf", width = 7, height = 7)
df.gapped.identical.dist.stats <- plotDensityOfInteractionRatios(df.gapped.identical.dist, 
                                                                 transform.fxn = log10, 
                                                                 add.legend = TRUE, 
                                                                 legend.labels = legend.labels$label, 
                                                                 legend.loc="topright"); 
df.gapped.identical.dist.stats 
dev.off()
## FROM COMPLETE UNGAPPED GRID
pdf("03b-complete-contactGrid.clusterSummitLog10FoldEnrichentOverFreqsOfIdenticalLengths.pdf", width = 7, height = 7)
df.complete.identical.dist.stats <- plotDensityOfInteractionRatios(df.complete.identical.dist[df.complete.identical.dist$FE <= 1e6,], 
                                                                   transform.fxn = log10, 
                                                                   add.legend = TRUE, 
                                                                   legend.labels = legend.labels$label, 
                                                                   legend.loc="topright"); df.complete.identical.dist.stats 
dev.off()

## PLOT RAW FOLD ENRICHMENT AND SHOW STATISTICS FROM BOTH (SHOULD BE THE SAME OR HIGHLY SIMILAR)
## FROM GAPPED GRID
pdf("04a-gapped-contactGrid.clusterSummitFoldEnrichentOverFreqsOfIdenticalLengths.FE_le_100.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.gapped.identical.dist[df.gapped.identical.dist$FE<=100,], 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, 
                               legend.labels = legend.labels$label, 
                               legend.loc="topright")
dev.off()
## FROM COMPLETE GRID
pdf("04b-complete-contactGrid.clusterSummitFoldEnrichentOverFreqsOfIdenticalLengths.FE_le_100.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.complete.identical.dist[df.complete.identical.dist$FE<=100, ], 
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
pdf("05-complete-contactGrid.clusterSummitLog10FoldEnrichentOverAllFreqsInGrid.pdf", width = 7, height = 7)
df.complete.entire.grid.stats <- plotDensityOfInteractionRatios(df.complete.entire.grid[df.complete.entire.grid$FE < 1e6, ], 
                                                                transform.fxn = log10, 
                                                                add.legend = TRUE, 
                                                                legend.labels = legend.labels$label, 
                                                                legend.loc="topleft"); df.complete.entire.grid.stats
dev.off()
## PLOT RAW FOLD ENRICHMENT 
pdf("06-complete-contactGrid.clusterSummitFoldEnrichentOverAllFreqsInGrid.FE_le_100.pdf", width = 7, height = 7)
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
# df.complete.short.range

pdf("07-complete-contactGrid.shortRangeIntFreqsFoldEnrichmentOverClusterSummitIntFreqs.pdf", width = 7, height = 7)
plotDensityOfInteractionRatios(df.complete.short.range, 
                               d.bw=5, transform.fxn = identity, xlab="Raw Fold Enrichment", 
                               add.legend = TRUE, legend.labels = legend.labels$label, legend.loc="topright")
dev.off()



## SCATTER PLOT OF ALL INTERACTIONS::  DISTANCES VS FREQUENCIES.

## LOG10 FREQS AND DISTS (LOG LOG PLOT) - SCALED BY SHORT RANGE
## USE ONLY FREQS > 0; i.e. NONZERO.
## FILTERING TINY VALUES FURTHER IMPROVES PLOT BY COLLAPSING Y-AXIS A BIT.
## TAKE THE TOP 99.95% OF NON-ZERO VALUES (EXCLUDE BOTTOM 0.05% OF INT FREQS > 0)
pdf("08a-complete-contactGrid.scatterPlot.IntDists-vs-IntFreqs.freqsScaledByMedShortRange.log10X.log10Y.downsample100k.pdf", width = 7, height = 7)
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

pdf("08b-complete-contactGrid.scatterPlot.IntDists-vs-IntFreqs.freqsScaledByMedShortRange.log10X.log10Y.all-data.pdf", width = 7, height = 7)
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


########################
## OTHER EXAMPLES

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

