##############################################################################################
### SCOPE - Scatter Cluster Of Paired Ends - Analysis Pipeline.
###       - Version v1.20250307 ; SCOPE-v1.20250307.R
### FUNCTIONS. Simply use source("/Path/To/This/File.R") in an R notebook.
###            - Recommended: start with one of our templates R scripts in Rstudio.
###                         : Aim it at your pre-processed text file (see below)
###                         : Then do cmd+enter through it.
##############################################################################################
### Cite the following pre-print (or follow-up publication associated with it):
###    Chromosome-scale scaffolding of the fungus gnat genome (Diptera: Bradysia coprophila)
###    John M. Urban, Susan A. Gerbi, Allan C. Spradling
###    bioRxiv 2022.11.03.515061; doi: https://doi.org/10.1101/2022.11.03.515061
##############################################################################################
### The first published version, used in the published paper, is part of the code section of:
###   - https://github.com/JohnUrban/Bcop_v2
### Updates to it will be hosted at: https://github.com/JohnUrban/SCOPE
##############################################################################################
### What you will need to run SCOPE:
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
### VERSION NOTES:
### - Current Version:  v1.20250307
### - Previous versions: v1.20250305, v0.20241023
##############################################################################################
### FUNCTIONS BELOW:
### - read.extra.short.pe
### - contactGrid.creation
### - contactGrid.cell.peak.calling
### - return.PeakCellInfo.From.Different.ContactGrid
### - contactGrid.clusterPeakCells.CallRegions
### - contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios
### - contactGrid.clusterPeakCells.RegionStats 
### - contactGrid.clusterPeakCells.CallRegions.hierarchical
### - contactGrid.clusterPeakCells.CallRegions.kmeans
### - contactGrid.ClusteredRegions.Filter
### - plotClusterDendrogram
### - return.ClusterInfo.From.Different.ContactGrid
### - min.max.normalize
### - log10.min.max.normalize
### - median.normalize
### - log10.median.normalize
### - z.score.normalize 
### - getContactGridLatticePlotParameters
### - plotContactGrid.Raw
### - plotContactGrid.WithAllCellPeaks
### - plotContactGrid.WithClusteredPeakRegions
### - vectorizeIntFreqsAndIntDists
### - compareClusterContactFreqSummitsToBackground
### - getShortRangeIntFreqs 
### - compareClusterContactFreqSummitsToShortRange
### - makeLabelsForClusters
### - plotDensityOfInteractionRatios
### - summarizeDF
### - contactGrid.scatterPlot














library(KernSmooth)
library(lattice)



## FUNCTIONS : READ EXTRASHORT.HIC.TXT FORMAT
read.extra.short.pe <- function(fh){
  read.table(fh, col.names = c("chr1","pos1","chr2", "pos2","mapq")) 
}


contactGrid.creation <- function(hic, tig_size, MAPQ=10, bx=250000, by=NA, gridsize=1024, minpair=1000e3, maxpair=NA){
  ## HIC is typically extra-short hic format, but actually it just needs to have two columns named "pos1" and "pos2".
  
  ## INTERPRET OPTS
  ## 1. Check bandwidth for y: "by"
  if(is.na(by)){by <- bx}
  
  ## FILTER XLINKS FOR MAPQ (SAVE IN HIC VAR). 
  ##  NOTE: IN TESTS, MAPQ FILTERING WORKED MUCH BETTER. 
  ##        NO FILTERING (MAPQ>=0) LEADS TO SEVERAL (OR MANY) EXTRA LIKELY FALSE PEAKS COMPARED TO MAPQ>=10.
  if(MAPQ>0){ hic <- hic[hic$mapq>=MAPQ, ] }
  
  ## PARE HIC OBJ DOWN TO TWO COLUMNS
  hic <- hic[,c("pos1","pos2")]
  
  ## DEFINE MAX DISTANCE BETWEEN MATES.
  ## Check maxpair, if NA, then make tig_size.
  #maxpair <- tig_size
  if(is.na(maxpair)){maxpair <- tig_size} ## this can't be longer than contig or chrom ; causes un-alerted error in plot.
  pairsizerange <- c(minpair, maxpair)
  
  ## DEFINE BOUNDS OF CONTIG (for now just using all)
  contig_range <- c(0, tig_size)
  
  ## KEEP ONLY MATE PAIRS WITH SPECIFIED DISTANCES BETWEEN MATES ACCORDING TO MINPAIR AND MAXPAIR.
  hicfrags <- hic[abs(hic$pos1-hic$pos2)>=minpair & abs(hic$pos1-hic$pos2)<=maxpair, ]
  
  ## CREATE MIRROR ACROSS DIAGONAL 
  revpairs <- hicfrags
  revpairs$pos1 <- hicfrags$pos2
  revpairs$pos2 <- hicfrags$pos1
  ## MIRROR NEEDS TO BE REVISITED: IS THIS "OVERCOUNTING" THE SELF-SELF DIAGONAL....?
  
  ## DEFINE ENTIRE SCATTER PLOT COMBINING MATES AND MIRRORED MATES; AND SORT THEM.
  xy <- rbind(hicfrags, revpairs)
  hicfrags <- xy[order(xy$pos1, xy$pos2),]
  
  ## GET DENSITY VALUES ACROSS SCATTER PLOT TO GET THE CONTACT GRID
  contactGrid <- bkde2D(x = hicfrags, 
                        bandwidth = c(bx,by), 
                        gridsize = c(gridsize, gridsize), 
                        range.x = list(x1=contig_range, 
                                       x2=contig_range))
  
  ## ADD INFO ON GRID CELL WIDTHS IN BP
  N <- length(contactGrid$x1)
  M <- N-1
  contactGrid$grid.cell.length <- median(as.vector(contactGrid$x1[2:N]) - contactGrid$x1[1:M] )
  
  ## RETURN THE CONTACT GRID
  return( contactGrid ) 
}

contactGrid.cell.peak.calling <- function(contactGrid, 
                                          quantile_thresh=0.998){
  
  ## GRID CELL PEAK CALLING.
  ## NOTE: x=contactGrid$x1; y=contactGrid$x2; z=contactGrid$fhat
  
  ## DEFINE PEAK THRESHOLD
  thresh <- quantile(contactGrid$fhat, 
                     quantile_thresh) #Good values: 0.999 , 0.9987*, 0.9985, 0.998
  
  ## ITERATE OVER CELLS IN GRID TO EXTRACT XY-COORDS OF CELLS THAT EXCEED THRESHOLD.
  xpeak<-c()
  ypeak<-c()
  zpeak <- c()
  grid.x <- c()
  grid.y <- c()
  for(i in 1:length(contactGrid$x1)){
    for(j in i:length(contactGrid$x2)){
      if(contactGrid$fhat[i,j]>thresh){
        xpeak  <- c(xpeak, contactGrid$x1[i])
        ypeak  <- c(ypeak, contactGrid$x2[j])
        zpeak  <- c(zpeak, contactGrid$fhat[i,j])
        grid.x <- c(grid.x, i)
        grid.y <- c(grid.y, j)
      }
    }
  }
  
  ## CREATE PEAK COORDS VAR FROM RESULTS
  peaks <- data.frame(x=xpeak,
                      y=ypeak,
                      z=zpeak,
                      grid.x=grid.x,
                      grid.y=grid.y)
  
  ## RETURN PEAK COORDS
  return(peaks)
}


return.PeakCellInfo.From.Different.ContactGrid <- function(peaks, contactGrid){
  new.peaks <- peaks
  for(i in 1:nrow(peaks)){
    grid.x <- peaks$grid.x[i]
    grid.y <- peaks$grid.y[i]
    new.peaks$x[i] <- contactGrid$x1[grid.x]
    new.peaks$y[i] <- contactGrid$x2[grid.y]
    new.peaks$z[i] <- contactGrid$fhat[grid.x, grid.y]
  }
  return(new.peaks)
} 




contactGrid.clusterPeakCells.CallRegions <- function(peaks,
                                                     k, 
                                                     method="h",
                                                     iter.max = 1000, nstart = 1000,
                                                     seglwd=0.5,
                                                     segcol="black"){
  ## METHODS: only looks at first letter. h for hierarchical, k for kmeans. for hierarchical clustering use any word that starts with "h"; for k-means, any that starts with k. 
  method <- substr(method,start=1, stop=1)
  if(method=="h"){
    contactGrid.clusterPeakCells.CallRegions.hierarchical(peaks=peaks,
                                                          k=k,
                                                          seglwd=seglwd,
                                                          segcol=segcol)
  }else if(method == "k"){
    contactGrid.clusterPeakCells.CallRegions.kmeans (peaks=peaks,
                                                     k=k,
                                                     iter.max = iter.max, nstart = nstart,
                                                     seglwd=seglwd,
                                                     segcol=segcol)
  }
  
}


contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios <- function(x){
  return(x/min(x))
}





contactGrid.clusterPeakCells.RegionStats <- function(peaks, groups, k=NA){
  if(is.na(k)){k <- max(groups)}
  hx <- c() 
  hy <- c()
  hxmin <- c()
  hxmax <- c()
  hymin <- c()
  hymax <- c()
  zSummit <- c()
  xSummit <- c()
  ySummit <- c()
  z.med <- c()
  z.mad <- c()
  z.mean <- c()
  z.stdev <- c()
  z.min <- c()
  z.max <- c()
  zSummit.grid.x <- c()
  zSummit.grid.y <- c()
  for (i in 1:k){
    hx <- c(hx, median(peaks$x[groups == i]))
    hy <- c(hy, median(peaks$y[groups == i]))
    #hx<- c(hx, mean(peaks$x[groups == i]))
    #hy<- c(hy, mean(peaks$y[groups == i]))
    ###########
    hxmin <- c(hxmin, min(peaks$x[groups == i]))
    hxmax <- c(hxmax, max(peaks$x[groups == i]))
    hymin <- c(hymin, min(peaks$y[groups == i]))
    hymax <- c(hymax, max(peaks$y[groups == i]))
    
    summitIndex <- which.max(peaks$z[groups == i])
    zSummit <- c(zSummit, peaks$z[groups == i][summitIndex])
    xSummit <- c(xSummit, peaks$x[groups == i][summitIndex])
    ySummit <- c(ySummit, peaks$y[groups == i][summitIndex])
    zSummit.grid.x <- c(zSummit.grid.x, peaks$grid.x[groups == i][summitIndex])
    zSummit.grid.y <- c(zSummit.grid.y, peaks$grid.y[groups == i][summitIndex])
    med <- median(peaks$z[groups == i])
    mad <- median( abs(peaks$z[groups == i] - med) )
    z.med <- c(z.med, med)
    z.mad <- c(z.mad, mad )
    z.mean <- c(z.mean, mean(peaks$z[groups == i]))
    z.stdev <- c(z.stdev, sd(peaks$z[groups == i]))
    z.min <- c(z.min, min(peaks$z[groups == i]))
    z.max <- c(z.max, max(peaks$z[groups == i]))

    
  }
  
  
  z.med.ratio <- contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios( z.med )
  z.mean.ratio <- contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios( z.mean )
  z.summit.ratio <- contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios( zSummit )
  z.max.ratio <- contactGrid.clusterPeakCells.RegionStats.helperFxn.Ratios( z.max )

  
  cluster.stats <- list(hx=hx,
                        hy=hy,
                        hxmin=hxmin,
                        hxmax=hxmax,
                        hymin=hymin,
                        hymax=hymax,
                        xSummit=xSummit,
                        ySummit=ySummit,
                        zSummit=zSummit,
                        z.med=z.med,
                        z.mad=z.mad,
                        z.mean=z.mean,
                        z.stdev=z.stdev,
                        z.min=z.min,
                        z.max=z.max,
                        z.med.ratio=z.med.ratio,
                        z.mean.ratio=z.mean.ratio,
                        z.summit.ratio=z.summit.ratio,
                        z.max.ratio=z.max.ratio,
                        zSummit.grid.x=zSummit.grid.x,
                        zSummit.grid.y=zSummit.grid.y)
  
  return(cluster.stats)
}



contactGrid.clusterPeakCells.CallRegions.hierarchical <- function(peaks,
                                                     k,
                                                     seglwd=0.5,
                                                     segcol="black"){
  ### HCLUST TO CALL PEAKS
  distmat <- dist(peaks) 
  hier <- hclust(distmat)
  groups <- cutree(hier, k = k)
  
  ## GET CLUSTER STATS
  cluster.stats <- contactGrid.clusterPeakCells.RegionStats(peaks = peaks, 
                                                            groups = groups,
                                                            k = k)
  
  ## DEFINE CLUSTERS LIST
  clusters <- list(distmat=distmat,
                   hier=hier,
                   groups=groups,
                   hx=cluster.stats$hx,
                   hy=cluster.stats$hy,
                   hxmin=cluster.stats$hxmin,
                   hxmax=cluster.stats$hxmax,
                   hymin=cluster.stats$hymin,
                   hymax=cluster.stats$hymax,
                   xSummit=cluster.stats$xSummit,
                   ySummit=cluster.stats$ySummit,
                   zSummit=cluster.stats$zSummit,
                   z.med=cluster.stats$z.med,
                   z.mad=cluster.stats$z.mad,
                   z.mean=cluster.stats$z.mean,
                   z.stdev=cluster.stats$z.stdev,
                   z.min=cluster.stats$z.min,
                   z.max=cluster.stats$z.max,
                   z.med.ratio=cluster.stats$z.med.ratio,
                   z.mean.ratio=cluster.stats$z.mean.ratio,
                   z.summit.ratio=cluster.stats$z.summit.ratio,
                   z.max.ratio=cluster.stats$z.max.ratio,
                   zSummit.grid.x=cluster.stats$zSummit.grid.x,
                   zSummit.grid.y=cluster.stats$zSummit.grid.y,
                   hx.center=NA,
                   hy.center=NA,
                   kmeans=NA)
  
  ## RETURN CLUSTERS
  return(clusters)
}




contactGrid.clusterPeakCells.CallRegions.kmeans <- function(peaks,
                                                            k, iter.max = 1000, nstart = 1000,
                                                            seglwd=0.5,
                                                            segcol="black"){
  ### KMEANS TO CALL PEAKS
  clustersK <- kmeans(peaks, centers=k, iter.max = iter.max, nstart = nstart)
  


  
  ## GET CLUSTER STATS
  ## To make data structure same as for hierarchical clust output.
  groups <- clustersK$cluster
  cluster.stats <- contactGrid.clusterPeakCells.RegionStats(peaks = peaks, 
                                                            groups = groups,
                                                            k = k)
  
  ## Extra info specific to kmeans
  hx.center <- as.vector(clustersK$centers[,1])
  hy.center <- as.vector(clustersK$centers[,2])

  
  
  ## DEFINE CLUSTERS LIST
  clusters <- list(distmat=NA,
                   hier=NA,
                   groups=groups,
                   hx=cluster.stats$hx,
                   hy=cluster.stats$hy,
                   hxmin=cluster.stats$hxmin,
                   hxmax=cluster.stats$hxmax,
                   hymin=cluster.stats$hymin,
                   hymax=cluster.stats$hymax,
                   xSummit=cluster.stats$xSummit,
                   ySummit=cluster.stats$ySummit,
                   zSummit=cluster.stats$zSummit,
                   z.med=cluster.stats$z.med,
                   z.mad=cluster.stats$z.mad,
                   z.mean=cluster.stats$z.mean,
                   z.stdev=cluster.stats$z.stdev,
                   z.min=cluster.stats$z.min,
                   z.max=cluster.stats$z.max,
                   z.med.ratio=cluster.stats$z.med.ratio,
                   z.mean.ratio=cluster.stats$z.mean.ratio,
                   z.summit.ratio=cluster.stats$z.summit.ratio,
                   z.max.ratio=cluster.stats$z.max.ratio,
                   zSummit.grid.x=cluster.stats$zSummit.grid.x,
                   zSummit.grid.y=cluster.stats$zSummit.grid.y,
                   hx.center=hx.center,
                   hy.center=hy.center,
                   kmeans=clustersK)
  
  ## RETURN CLUSTERS
  return(clusters)
}



contactGrid.ClusteredRegions.Filter <- function(clusters, minRegionSeparation=5e6, updateClusters=FALSE){

  if(updateClusters){
    gate <- abs( clusters$hx-clusters$hy ) > minRegionSeparation
    clusters$hx <- clusters$hx[gate]
    clusters$hy <- clusters$hy[gate]
    clusters$hxmin <- clusters$hxmin[gate]
    clusters$hxmax <- clusters$hxmax[gate]
    clusters$hymin <- clusters$hymin[gate]
    clusters$hymax <- clusters$hymax[gate]
    clusters$xSummit <- clusters$xSummit[gate]
    clusters$ySummit <- clusters$ySummit[gate]
    clusters$zSummit <- clusters$zSummit[gate]
    clusters$z.med <- clusters$z.med[gate]
    clusters$z.mad <- clusters$z.mad[gate]
    clusters$z.mean <- clusters$z.mean[gate]
    clusters$z.stdev <- clusters$z.stdev[gate]
    clusters$z.min <- clusters$z.min[gate]
    clusters$z.max <- clusters$z.max[gate]
    clusters$z.med.ratio <- clusters$z.med.ratio[gate]
    clusters$z.mean.ratio <- clusters$z.mean.ratio[gate]
    clusters$z.summit.ratio <- clusters$z.summit.ratio[gate]
    clusters$z.max.ratio <- clusters$z.max.ratio[gate]
    clusters$zSummit.grid.x <- clusters$zSummit.grid.x[gate]
    clusters$zSummit.grid.y <- clusters$zSummit.grid.y[gate]
    return(clusters)
  }else{
    ust <- clusters$hx
    uen <- clusters$hy
    ust <- ust[abs( clusters$hx-clusters$hy ) > minRegionSeparation] #; ust
    uen <- uen[abs( clusters$hx-clusters$hy ) > minRegionSeparation] #; uen
    return(data.frame(ust=ust, 
                      uen=uen))
  }
}




plotClusterDendrogram <- function(cluster, method="h"){
  if(method!="h"){
    print("Only for hieraechical clustering at the moment.")
  }else{
    plot(clusters$hier)
  }
}



return.ClusterInfo.From.Different.ContactGrid <- function(peaks, contactGrid, clusters){
  new.peaks <- return.PeakCellInfo.From.Different.ContactGrid(peaks, contactGrid)
  cluster.stats <- contactGrid.clusterPeakCells.RegionStats(new.peaks, clusters$groups)
  return(cluster.stats)
} 


min.max.normalize <- function(z){
  return(
    (z-min(z))/(max(z)-min(z))
  )
}



log10.min.max.normalize <- function(z){
  return(
    log10( min.max.normalize(z) )
  )
}

median.normalize <- function(z, add.nonzero.min.pseudo=FALSE){
  if(add.nonzero.min.pseudo){ z <- z + min( z[z>0] ) }
  return(
    z/median(z)
  )
}

log10.median.normalize <- function(z){
  ## Assumes all values > 0.
  return(
    log10( median.normalize(z, add.nonzero.min.pseudo=TRUE  ) )
  )
}


z.score.normalize <- function(z){
  return(
    (z-mean(z))/sd(z)
  )
}



getContactGridLatticePlotParameters <- function(contactGrid, defaultminz=0, 
                                                coltrio=c("dark blue","white","red"), 
                                                transform.fxn=identity, sigdig=3){
  zvals <- transform.fxn( contactGrid$fhat )
  # Get maximum value in Z dimension
  maxz <- max( zvals )
  minz <- min( zvals )
  zwidth <- maxz - minz
  
  # compute these outside of list as they are needed by colorkey
  # seqmaxz  = seq( min( c( defaultminz, minz ) ), 
  #                 maxz, 
  #                 abs(maxz/1000))
  # seqmaxz2 = seq( min( c( defaultminz, minz ) ), 
  #                 maxz, 
  #                 abs(maxz/3))
  seqmaxz  = seq( minz, 
                  maxz, 
                  abs(zwidth/1000))
  seqmaxz2 = seq( minz, 
                  maxz, 
                  abs(zwidth/4))
  print(c("Z-range:",minz, maxz, zwidth))
  ## Fill up list with parameters needed.
  contactGridLatticePlotParameters <- list(x        = contactGrid$x1,
                                           y        = contactGrid$x2,
                                           z        = zvals,
                                           maxz     = maxz,
                                           seqmaxz  = seqmaxz,
                                           seqmaxz2 = seqmaxz2,
                                           colorkey = list(space = "right", 
                                                           col   = colorRampPalette(coltrio), 
                                                           at    = seqmaxz, 
                                                           raster= TRUE, 
                                                           tck   = 1, 
                                                           labels= list(labels = round( seqmaxz2, digits = sigdig ), 
                                                                        at     = seqmaxz2 )))
  # Return list.
  return(contactGridLatticePlotParameters)
}




plotContactGrid.Raw <- function(contactGrid, seglwd=0.5, segcol="black", coltrio=c("dark blue","white","red"), 
                                transform.fxn=identity, defaultminz=0, sigdig=3){
  
  ## get Contact Grid Lattice Plot Parameters
  contactGridLatticePlotParameters <- getContactGridLatticePlotParameters( contactGrid=contactGrid, 
                                                                           defaultminz=defaultminz, 
                                                                           coltrio=coltrio, 
                                                                           transform.fxn=transform.fxn,
                                                                           sigdig=sigdig )
  
  # ## WITH SQUARES ALONG DIAGONAL
  # ust <- peaks$x
  # uen <- peaks$y
  
  ## PLOT
  # lattice.options(axis.padding=list(factor=0.05))
  
  levelplot(contactGridLatticePlotParameters$z, 
            row.values    = contactGridLatticePlotParameters$x, 
            column.values = contactGridLatticePlotParameters$y, 
            colorkey      = contactGridLatticePlotParameters$colorkey, 
            xlab="", ylab="", 
            scales=list(x=list(rot=90), 
                        pretty=TRUE, ylab="", 
                        xlab="", 
                        tck = c(1,0)),
            col.regions = colorRampPalette(coltrio), 
            at = contactGridLatticePlotParameters$seqmaxz)
}



plotContactGrid.WithAllCellPeaks <- function(contactGrid, peaks, seglwd=0.5, segcol="black", 
                                             coltrio=c("dark blue","white","red"), 
                                             transform.fxn=identity, defaultminz=0, sigdig=3){
  
  ## get Contact Grid Lattice Plot Parameters
  contactGridLatticePlotParameters <- getContactGridLatticePlotParameters( contactGrid=contactGrid, 
                                                                           defaultminz=defaultminz, 
                                                                           coltrio=coltrio, 
                                                                           transform.fxn=transform.fxn,
                                                                           sigdig=sigdig )
  
  ## WITH SQUARES ALONG DIAGONAL
  ust <- peaks$x
  uen <- peaks$y
  
  ## PLOT
  # lattice.options(axis.padding=list(factor=0.05)) ## Probably need this.
  levelplot(contactGridLatticePlotParameters$z, 
            row.values    = contactGridLatticePlotParameters$x, 
            column.values = contactGridLatticePlotParameters$y, 
            colorkey      = contactGridLatticePlotParameters$colorkey, 
            xlab="", ylab="", 
            scales = list(x=list(rot=90), 
                          pretty=TRUE, ylab="", 
                          xlab="", 
                          tck = c(1,0)),
            col.regions = colorRampPalette(coltrio), 
            at = contactGridLatticePlotParameters$seqmaxz,
            panel = function(...){panel.levelplot.raster(...); 
              panel.segments(x0=c(ust, ust, uen, uen), 
                             x1=c(ust, uen, uen, ust), 
                             y0=c(ust, ust, uen, uen), 
                             y1=c(uen, ust, ust, uen), 
                             lwd=seglwd, 
                             col=segcol)})
}



plotContactGrid.WithClusteredPeakRegions <- function(contactGrid, clusters, minRegionSeparation=5e6,
                                                     seglwd=0.5, segcol="black", 
                                                     coltrio=c("dark blue","white","red"),
                                                     transform.fxn=identity, defaultminz=0, sigdig=3){
  
  ## get Contact Grid Lattice Plot Parameters
  contactGridLatticePlotParameters <- getContactGridLatticePlotParameters( contactGrid=contactGrid, 
                                                                           defaultminz=defaultminz, 
                                                                           coltrio=coltrio, 
                                                                           transform.fxn=transform.fxn,
                                                                           sigdig=sigdig )
  
  ## WITH SQUARES ALONG DIAGONAL
  squares <- contactGrid.ClusteredRegions.Filter(clusters=clusters, 
                                                 minRegionSeparation=minRegionSeparation)
  ust <- squares$ust
  uen <- squares$uen
  
  ## PLOT
  levelplot(contactGridLatticePlotParameters$z, 
            row.values    = contactGridLatticePlotParameters$x, 
            column.values = contactGridLatticePlotParameters$y, 
            colorkey      = contactGridLatticePlotParameters$colorkey, 
            xlab="", ylab="", 
            scales=list(x=list(rot=90), 
                        pretty=TRUE, ylab="", 
                        xlab="", 
                        tck = c(1,0)),
            col.regions = colorRampPalette(coltrio), 
            at = contactGridLatticePlotParameters$seqmaxz,
            panel = function(...){panel.levelplot.raster(...); 
              panel.segments(x0=c(ust, ust, uen, uen), 
                             x1=c(ust, uen, uen, ust), 
                             y0=c(ust, ust, uen, uen), 
                             y1=c(uen, ust, ust, uen), 
                             lwd = seglwd, 
                             col = segcol)})
}








#######################################################################################


vectorizeIntFreqsAndIntDists <- function(contactGrid){
  #' vectorizeIntFreqsAndIntDists
  #' Takes in contactGrid (an  output of contactGrid.creation() or bkde2D())
  #' - Assumes contactGrid is square: ncols = nrows; i.e. NxN.
  #' Updates it with two new related "columns" (as if a data.frame)
  #' 1. IntDists : Interaction Distance for the given Interaction Frequency with the same index.
  #' 2. IntFreqs : Interaction Frequency for the given Interaction Distance with the same index.

  
  ## Extract grid dimensions
  gridsize <- length(contactGrid$x1)
  nrows <- gridsize
  ncols <- gridsize
  
  ## Initialize vectors to store variables.
  nUniqueXYcombinations <- sum(1:gridsize)
  contactGrid$IntDists <- numeric(length = nUniqueXYcombinations) ## actual vector length will be sum(1:gridsize) i.e. nUniqueXYcombinations
  contactGrid$IntFreqs <- numeric(length = nUniqueXYcombinations)  ## actual vector length will be sum(1:gridsize) i.e. nUniqueXYcombinations
  
  ## Initialize index variable, idx
  idx = 0
  for( i in 1:nrows){
    for(j in i:ncols){
      idx <- idx+1
      contactGrid$IntDists[idx] <- abs( contactGrid$x1[i] - contactGrid$x2[j] )
      contactGrid$IntFreqs[idx] <- contactGrid$fhat[i,j]
    }
  }
  return(contactGrid)
}



compareClusterContactFreqSummitsToBackground <- function(clusters, contactGrid, distBufferFactor=0, scaleByBufferFactor=FALSE, enforce.nonzero.intfreq=FALSE){
  #' Takes in:
  #' - clusters, 
  #' - contactGrid, (an  output of contactGrid.creation() or bkde2D())
  #'    - Should be updated with vectorizeIntFreqsAndIntDists()
  #'    - But, this will check and run it if need be.
  #' - distBufferFactor, the max difference from the interaction distance of interest.
  #'    - Used to define allowable interaction distances to interrogate for interaction frequencies.
  #'    - Default = 0.
  #'    -   i.e. only identical interaction distances to the one of interest.
  #'    -   Some times there are not enough values when requiring the same value.
  #'    -   Thus using similar distances, as opposed to identical, is useful to sample more data.
  #'    - Setting this to 100e3, for example, would allow distances to be the IntDist +/- 100kb.
  #'    - Setting this to something very high, such as 100000000000000, compares it to ALL INTERACTION DISTANCES.
  #' - scaleByBufferFactor, whether to use distBufferFactor as a max difference or to scale the current IntDist to get the max difference.
  #'    -   By default (FALSE), all distances (X) with an absolute difference of X-currIntDist <= distBufferFactor will be used.
  #'    -   If this is set (TRUE), then _ will use all distances with an abs diff <= distBufferFactor*currIntDist
  #'        - When this option used, values are typically between 0-1.
  #'        - i.e. if you want to collect all distances that are +/- 5% or 10%, then set distBufferFactor to  0.05 or 0.1 and set this as TRUE.
  #' Outputs:
  #'  - 2-column data.frame that can be used for plotting or otherwise.
  #'      - column 1 is the cluster index.
  #'      - column 2 is the cluster summit FE values.
  #' Terminology found within:
  #'    1. IntDists : Vector of Interaction Distances   for Interaction Frequencies with the same index.
  #'    2. IntFreqs : Vector of Interaction Frequencies for Interaction Distances   with the same index.
  
  ## Check contactGrid for IntFreqs and IntDists; if not present, update with them.
  if ( ! 'IntDists' %in% names(contactGrid) ){
    contactGrid <- vectorizeIntFreqsAndIntDists(contactGrid)
  }
  
  ## Extract interaction frequency point estimates (IntFreqs) from clusters (i.e. "Dots" on Hi-C maps).
  dot.IntFreq <- clusters$zSummit
  
  ## Extract interaction distances (IntDists) from clusters.
  dot.IntDist <- clusters$ySummit - clusters$xSummit
  
  ## Compute Fold-Enrichment (FE) of cluster IntFreqs to background IntFreqs of similar IntDists.
  ## Initializw
  all.dots.idx <- c()
  all.dots.FE <- c()
  nClusters <- length(clusters$zSummit)
  ## For each cluster, get IntFreq_Clust / IntFreq_x (for x in all IntFreqs of similar IntDists)
  for(i in 1:nClusters){
    ## Get BOOL INDEX of 
    if(scaleByBufferFactor){
      IntDistGate <- abs( contactGrid$IntDists - dot.IntDist[i] ) <= distBufferFactor * dot.IntDist[i]
    } else {
      IntDistGate <- abs( contactGrid$IntDists - dot.IntDist[i] ) <= distBufferFactor
    }
    if(enforce.nonzero.intfreq){
      IntDistGate <- IntDistGate & contactGrid$IntFreqs > 0
    }
    dot.intFreq.FE <- dot.IntFreq[i] / contactGrid$IntFreqs[ IntDistGate ]
    all.dots.FE  <- c(all.dots.FE, dot.intFreq.FE)
    all.dots.idx <- c(all.dots.idx, rep(i, length(dot.intFreq.FE)))

  }
  all.dots.FE.df <- data.frame(idx=all.dots.idx, FE=all.dots.FE)
  return(all.dots.FE.df)
}




getShortRangeIntFreqs <- function(contactGrid, maxGridCellDist=1, interpretAsNtDistance=FALSE, enforce.nonzero.intfreq=FALSE){
  ## Get BOOL INDEX of 
  if(interpretAsNtDistance){
    IntDistGate <- abs( contactGrid$IntDists ) <= maxGridCellDist
  } else {
    IntDistGate <- abs( contactGrid$IntDists ) <= maxGridCellDist * contactGrid$grid.cell.length
  }
  if(enforce.nonzero.intfreq){
    IntDistGate <- IntDistGate & contactGrid$IntFreqs > 0
  }
  
  ## Extract short-range distances.
  shortRangeIntFreqs <- contactGrid$IntFreqs[ IntDistGate ]
  
  return(shortRangeIntFreqs)
}



compareClusterContactFreqSummitsToShortRange <- function(clusters, contactGrid, 
                                                         maxGridCellDist=1, interpretAsNtDistance=FALSE, 
                                                         enforce.nonzero.intfreq=FALSE){
  #' Takes in:
  #' - clusters, 
  #' - contactGrid, (an  output of contactGrid.creation() or bkde2D())
  #'    - Should be updated with vectorizeIntFreqsAndIntDists()
  #'    - But, this will check and run it if need be.
  #' - maxGridCellDist, the max interaction distance in units of grid cells to include in the "short range" group.
  #'    - Used to define allowable interaction distances to interrogate for interaction frequencies.
  #'    - Default = 1.
  #'    -   Use the number of interaction frequencies within the distance of 1 grid cell. i.e. only the diagonal.
  #' - interpretAsNtDistance, 
  #'    -   By default (FALSE), maxGridCellDist will be a multiple of contactGrid$grid.cell.length.
  #'    -   If this is set (TRUE), then will just be interpreted as the nucleotide (Nt) length.
  #'        - Example 1: maxGridCellDist=2, interpretAsNtDistance=FALSE; short-range will be all <= 2*contactGrid$grid.cell.length
  #'        - Example 2: maxGridCellDist=100e3, interpretAsNtDistance=TRUE; short-range will be all <= 100kb.
  #' Outputs:
  #'  - 2-column data.frame that can be used for plotting or otherwise.
  #'      - column 1 is the cluster index.
  #'      - column 2 is the cluster summit FE values.
  #' Terminology found within:
  #'    1. IntDists : Vector of Interaction Distances   for Interaction Frequencies with the same index.
  #'    2. IntFreqs : Vector of Interaction Frequencies for Interaction Distances   with the same index.
  
  ## Check contactGrid for IntFreqs and IntDists; if not present, update with them.
  if ( ! 'IntDists' %in% names(contactGrid) ){
    contactGrid <- vectorizeIntFreqsAndIntDists(contactGrid)
  }
  
  ## Extract interaction frequency point estimates (IntFreqs) from clusters (i.e. "Dots" on Hi-C maps).
  dot.IntFreq <- clusters$zSummit
  
  ## Extract interaction distances (IntDists) from clusters.
  dot.IntDist <- clusters$ySummit - clusters$xSummit
  
  ## Compute Fold-Enrichment (FE) of cluster IntFreqs to background IntFreqs of similar IntDists.
  ## Initializw
  all.dots.idx <- c()
  all.dots.FE <- c()
  nClusters <- length(clusters$zSummit)
  
  ## Get BOOL INDEX of 
  # if(interpretAsNtDistance){
  #   IntDistGate <- abs( contactGrid$IntDists ) <= maxGridCellDist
  # } else {
  #   IntDistGate <- abs( contactGrid$IntDists ) <= maxGridCellDist * contactGrid$grid.cell.length
  # }
  # if(enforce.nonzero.intfreq){
  #   IntDistGate <- IntDistGate & contactGrid$IntFreqs > 0
  # }
  ## Extract short-range distances once (instead of in every loop)
  #shortRangeIntFreqs <- contactGrid$IntFreqs[ IntDistGate ]
  
  ## Extract short-range distances once (instead of in every loop)
  shortRangeIntFreqs <- getShortRangeIntFreqs(contactGrid=contactGrid, 
                                              maxGridCellDist=maxGridCellDist, 
                                              interpretAsNtDistance=interpretAsNtDistance, 
                                              enforce.nonzero.intfreq=enforce.nonzero.intfreq)
  
  
  
  ## For each cluster, get IntFreq_Clust / IntFreq_x (for x in all IntFreqs of similar IntDists)
  for(i in 1:nClusters){
    dot.intFreq.FE <- shortRangeIntFreqs / dot.IntFreq[i]
    all.dots.FE  <- c(all.dots.FE, dot.intFreq.FE)
    all.dots.idx <- c(all.dots.idx, rep(i, length(dot.intFreq.FE)))
    
  }
  
  # Store as DF
  all.dots.FE.df <- data.frame(idx=all.dots.idx, FE=all.dots.FE)
  
  ## Return DF
  return(all.dots.FE.df)
}









makeLabelsForClusters <- function(clusters){
  nClusters <- length(clusters$zSummit)
  df <- data.frame(idx=1:nClusters, label=1:nClusters)
  for(i in 1:nClusters){
    df$label[i] <- paste(round(clusters$xSummit[i]), round(clusters$ySummit[i]), sep = "," )
  }
  return(df)
}




plotDensityOfInteractionRatios <- function(df, col1="idx", col2="FE", 
                                           d.from=NA, d.to=maxval, d.bw=0.05, d.n=1000, maxd.scale=0.05, 
                                           plot.title="", xlab="log10(Fold Enrichment)", ylab="Density",
                                           line.cols=c("black","red","blue","green","orange","purple","grey","gold","skyblue","magenta","yellow","violet","darkblue","darkgreen","cyan","darkcyan"),
                                           add.legend=FALSE, legend.labels=NA, legend.loc="topleft", legend.label.append.stat=TRUE,
                                           transform.fxn=log10, sigdig=4,
                                           vlines=NA){
  #' Takes in:
  #' - df, the "all.dots.FE.df" output from compareClusterContactFreqSummitsToBackground()
  #'    - or by default, any df that has "idx" and "FE" columns.
  #'    - change col1 and col2 variables for other columns.

  maxval <- max( transform.fxn( df[[col2]] ) )
  if(sum(is.na(d.from))>0){ 
    d.from <- min(0,  transform.fxn( df[[col2]] ) )  # if it is something tiny between 0-1, I just want it 0, but if its negative then use that.
  }
  density.list <- list()
  N <- max(df[[col1]])
  maxd.vec <- numeric(length=N)
  mind.vec <- numeric(length=N)
  for(i in 1:N){
    density.list[[i]] <-  density( transform.fxn( df[[col2]][df[[col1]]==i] ), from=d.from, to=maxval+d.bw, bw=d.bw, n=d.n)
    maxd.vec[i] <- max(density.list[[i]]$y)
    mind.vec[i] <- min(density.list[[i]]$y)
  }
  
  maxd <- max(maxd.vec)
  mind <- min(0, mind.vec) # if it is something tiny between 0-1, I just want it 0, but if its negative then use that.
  #print(c(d.from,maxval,mind,maxd))
  plot( density.list[[1]], ylim=c(mind, maxd+maxd.scale*maxd ), type="n" , xlab=xlab, ylab=ylab, main=plot.title, las=1)
  for(i in 1:N){
    lines(density.list[[i]], col=line.cols[i])
  }
  
  if(sum(is.na(vlines))==0){
    abline(v=vlines,lty=3)
  }
  
  statdf <- summarizeDF(df, col1, col2, sigdig, transform.fxn)
  
  if(add.legend){
    if(sum(is.na(legend.labels))>1){
      legend.labels <- 1:N
    }
    if(legend.label.append.stat){
      for(i in 1:N){
        legend.labels[i] <- paste(legend.labels[i],statdf$Med[i],sep=",")
      }
    }
    legend(legend.loc, legend=legend.labels, fill=line.cols[1:N])
  }
  
  return(statdf)
}



summarizeDF <- function(df, col1="idx", col2="FE", sigdig=2, transform.fxn=identity){
  N <- max(df[[col1]])
  statdf <- data.frame(Min=1:N, Q1=1:N, Med=1:N, Mean=1:N, Q3=1:N, Max=1:N)
  for(i in 1:N){
    statdf[i,] <- summary( round( transform.fxn( df[[col2]][df[[col1]]==i] ), sigdig ))
  }
  return(statdf)
}


contactGrid.scatterPlot <- function(contactGrid, clusters, tooShort=-1, 
                                    transform.fxn.freqs=identity, transform.fxn.dists=identity,
                                    down.sample=NA, scaleFreqsByMedShortRange=FALSE, grid.col="black", clust.col="red",
                                    maxGridCellDist=1, interpretAsNtDistance=FALSE, 
                                    enforce.nonzero.intfreq=FALSE, enforce.min.intfreq.quantile=FALSE, min.intfreq.quantile=0,
                                    grid.cex=0.2, clust.cex=1, grid.pch=19, clust.pch=1){
  #' tooShort, minimum Interaction Distance - 1.
  #'  - Will use all distances > tooShort.
  #'  - Default of -1 means use all distances >= 0 (i.e. >-1).
  #'  - Recommended to use something higher than 0 if the shortest distances are too high freq.
  #'  - Or use log10 freqs in that case
  #'  log10Dists and log10Freqs
  #'    - Use to log10 transform them before plotting. Can do one or both.
  #'  transform.fxn.freqs, transformation to apply to Interaction Frequencies.
  #'    - Typical values: 
  #'      - identity (default), 
  #'      - log10, 
  #'      - min.max.normalize
  #'      - median.normalize
  #'      - z.score.normalize
  #'      - log10.median.normalize
  #'      - log10.min.max.normalize
  #'      - Any function that processes a vector: e.g. sqrt
  #'  transform.fxn.freqs, transformation to apply to Interaction Frequencies.
  #'    - Typical values: 
  #'      - identity (default), 
  #'      - log10
  #'      - Other transformations for freqs can also be done, but not recommended.
  #' down.sample, use less data.
  #'  - default is to use all data points, but this can be massive and not plot well (millions of points)
  #'  - Provide an integer to down.sample -- sampling without replacement.
  #'  - If the integer is higher than total amount, then it will default to just using all; not sample with replacement.
 
  
  ## Check contactGrid for IntFreqs and IntDists; if not present, update with them.
  if ( ! 'IntDists' %in% names(contactGrid) ){
    contactGrid <- vectorizeIntFreqsAndIntDists(contactGrid)
  }
  
  ## INITIALIZE DIST AND FREQ VARIABLES
  Grid.IntFreqs <- contactGrid$IntFreqs
  Grid.IntDists <- contactGrid$IntDists
  Clust.IntFreqs <- clusters$zSummit
  Clust.IntDists <- clusters$ySummit - clusters$xSummit
  
  ## KEEP ONLY FREQS AND DISTS WHEN DIST > tooShort
  LengthBool <- Grid.IntDists > tooShort
  Grid.IntFreqs <- Grid.IntFreqs[LengthBool]
  Grid.IntDists <- Grid.IntDists[LengthBool]
  
  ## OPTIONAL FILTER: KEEP ONLY FREQS AND DISTS WHEN FREQ>0.
  if(enforce.nonzero.intfreq){
    FreqBool <- Grid.IntFreqs > 0
    Grid.IntFreqs <- Grid.IntFreqs[FreqBool]
    Grid.IntDists <- Grid.IntDists[FreqBool]
  }
  
  ## OPTIONAL FILTER: KEEP ONLY FREQS AND DISTS WHEN FREQ > min.quantile value.
  if(enforce.min.intfreq.quantile){
    min.intfreq <- as.vector(quantile(Grid.IntFreqs, min.intfreq.quantile))
    FreqBool <- Grid.IntFreqs > min.intfreq[1]
    Grid.IntFreqs <- Grid.IntFreqs[FreqBool]
    Grid.IntDists <- Grid.IntDists[FreqBool]
  }
  
  ## DOWN SAMPLE
  N <- length( Grid.IntDists )
  if(!is.na(down.sample)){
    if( down.sample > N ){
      down.sample <- N
    }
    indices <- sort(sample(x = 1:N, 
                           size = down.sample,
                           replace = FALSE))
  }else{
    indices <- 1:N
  }
  
  
  ## EXTRACT INT FREQS AND DISTS FROM GRID
  Grid.IntFreqs <- Grid.IntFreqs[indices]
  Grid.IntDists <- Grid.IntDists[indices]
  
  
  ## OPTIONAL NORMALIZE FREQS TO SHORT RANGE MEDIAN BEFORE PROCEEDING
  if(scaleFreqsByMedShortRange){
    shortRangeIntFreqs <- getShortRangeIntFreqs(contactGrid=contactGrid, 
                                                maxGridCellDist=maxGridCellDist, 
                                                interpretAsNtDistance=interpretAsNtDistance, 
                                                enforce.nonzero.intfreq=enforce.nonzero.intfreq) 
    medShortRangeVal <- median(shortRangeIntFreqs)
    Grid.IntFreqs  <- Grid.IntFreqs  / medShortRangeVal 
    Clust.IntFreqs <- Clust.IntFreqs / medShortRangeVal
  }
  
  ## FINAL TRANSFORMATIONS
  Grid.IntDists  <- transform.fxn.dists(Grid.IntDists)
  Grid.IntFreqs  <- transform.fxn.freqs(Grid.IntFreqs)
  Clust.IntDists <- transform.fxn.dists(Clust.IntDists)
  Clust.IntFreqs <- transform.fxn.freqs(Clust.IntFreqs)
  maxFreq <- max(Grid.IntFreqs, Clust.IntFreqs)
  minFreq <- min(0, Grid.IntFreqs, Clust.IntFreqs)
  
  print(c(minFreq,maxFreq))
  ## SCATTER PLOTTING
  plot(Grid.IntDists, 
       Grid.IntFreqs, 
       cex = grid.cex, 
       ylim=c(minFreq, maxFreq ), 
       col = grid.col,
       las = 1,
       pch=grid.pch)
  
  points(Clust.IntDists, 
         Clust.IntFreqs, 
         col=clust.col,
         cex=clust.cex,
         pch=clust.pch)
}
