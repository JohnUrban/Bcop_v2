## Set working directory
setwd("/Users/jurban/software/Bcop_v2/dev/EGGS/EGGS/")


## Download LAVE to ~/software : https://github.com/JohnUrban/lave
## Functions from the main Lave dir may change over time.
## To ensure this analysis is reproducible, a minimal set of code/functions from LAVE is also copied inside this directory (Bcop_v2/code/lave-dotplots)
## - The Lave code snippets here are fossilized as is, and will not undergo further dev.
source("~/software/lave/lave.R") 
source("~/software/lave//lave-dev.R") 
source("~/software/lave//lave-assembalign-dev.R") 



## System checks
system("pwd")
system("ls")
system("mkdir -p Rfigures")



### READING PAF FILES GENERATED IN "ncbi-dataset-processing" DIRECTORY OF Bcop_v2/code/EGGS (after steps 1-5)
cnames <- c("query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", 
            "match","alnlen","mapq","og","geneA","geneB","protA","protB","strandA","strandB")

bcop_dmel_paf <- read.table("../ncbi-dataset-processing/bcop_dmel.paf", 
                            as.is = TRUE, 
                            col.names = cnames)

bcop_aedes_paf <- read.table("../ncbi-dataset-processing/bcop_aedes.paf", 
                             as.is = TRUE, 
                             col.names = cnames)

bcop_anoph_paf <- read.table("../ncbi-dataset-processing/bcop_anoph.paf", 
                             as.is = TRUE, 
                             col.names = cnames)


bcop_bhyg_paf <- read.table("../ncbi-dataset-processing/bcop_bhyg-newnames.paf", 
                            as.is = TRUE, 
                            col.names = cnames)

bcop_bcop_paf <- read.table("../ncbi-dataset-processing/bcop_bcop.paf", 
                            as.is = TRUE, 
                            col.names = cnames[1:12])



## LOOKING AT PAR MAR
par(no.readonly=TRUE)$mar

par(mar=c(5.1,6.1,3.1,2.1)) ## This was fine for plotting in R, but needed changes for PDF (see below)



head(bcop_dmel_paf)

head(bcop_bhyg_paf)


## SOME PARAMETER SETTING
tgaps <- NA ;
qgaps <- NA ;
ngoncol=rgb(1,0,0,0.1); 
bgoncol=NA; 
qtigcol = rgb(0,0,0,0.35);
ttigcol = rgb(0,0,0,0.65)
gapcol<-"white"
circsize = 2.5; goncol=rgb(0,0,1,0.4); 




########################################################################################
########################################################################################
## DOT PLOTS BELOW
########################################################################################
########################################################################################



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PLOT DROSOPHILA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
paf <- bcop_dmel_paf
main<- "Drosophila melanogaster to Bradysia coprophila"; 
ylabels <- c("Bcop", "Dmel")
ttigcol = rgb(0.267004, 0.004874, 0.329415, 1);
qtigcol = rgb(0.647257, 0.8584, 0.209861, 1.0)
chrnames <- c("X","II","III","IV", "chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4") ## excluding chr Y
tnames = unique(paf[order(paf$tlen, decreasing = TRUE),6:7])$target
tassocnames = tnames[!(tnames %in% c("X","II","III","IV"))]
torder = c("X","II","III","IV", tassocnames)
qorder = c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrY", "chrMT", "NW_007931083.1")
qchr <- paf$query %in% chrnames
tchr <- paf$target %in% chrnames
chrpaf <- paf[qchr & tchr,]
dmelchrpaf <- chrpaf


pdf("Rfigures/bcop_dmel_dot-thicker-gridlines.pdf", width = 7, height = 7)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=dmelchrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", grid.subset.q = c(1,3,5), grid.lty = 1, qstep=20e6,
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, 
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
dev.off()

png("Rfigures/bcop_dmel_dot.png", width = 1000, height = 1000)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=dmelchrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", grid.subset.q = c(1,3,5),
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, 
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
dev.off()


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PLOT AEDES
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


paf <- bcop_aedes_paf
main<- "Aedes aegypti to Bradysia coprophila"; 
ylabels <- c("Bcop", "Aedes")
ttigcol = rgb(0.267004, 0.004874, 0.329415, 1);
qtigcol = rgb(0.185783, 0.704891, 0.485273, 1.0)
chrnames <- c("X","II","III","IV", "chr1", "chr2", "chr3")
tnames = unique(paf[order(paf$tlen, decreasing = TRUE),6:7])$target
tassocnames = tnames[!(tnames %in% c("X","II","III","IV"))]
torder = c("X","II","III","IV", tassocnames)
qorder = c("chr1", "chr2", "chr3", unique(paf$query)[4:length(unique(paf$query))])
qchr <- paf$query %in% chrnames
tchr <- paf$target %in% chrnames
chrpaf <- paf[qchr & tchr,]
aedeschrpaf <- chrpaf


pdf("Rfigures/bcop_aedes_dot.pdf", width = 7, height = 7)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=aedeschrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", qstep=200e6, #grid.subset.q = c(1,3,5),
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, ,
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
dev.off()


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PLOT ANOPHELES
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


paf <- bcop_anoph_paf
main<- "Anopheles gambiae to Bradysia coprophila"; 
ylabels <- c("Bcop", "Anopheles")
ttigcol = rgb(0.267004, 0.004874, 0.329415, 1);
qtigcol = rgb(0.35236, 0.783011, 0.392636, 1.0)
chrnames <- c("X","II","III","IV", "chrX", "chr2L", "chr2R", "chr3L", "chr3R")
tnames = unique(paf[order(paf$tlen, decreasing = TRUE),6:7])$target
tassocnames = tnames[!(tnames %in% c("X","II","III","IV"))]
torder = c("X","II","III","IV", tassocnames)
qorder = c("chrX", "chr2L", "chr2R", "chr3L", "chr3R", unique(paf$query)[6:length(unique(paf$query))])
qchr <- paf$query %in% chrnames
tchr <- paf$target %in% chrnames
chrpaf <- paf[qchr & tchr,]
anophchrpaf <- chrpaf

pdf("Rfigures/bcop_anopheles_dot.pdf", width = 7, height = 7)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=anophchrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", grid.subset.q = c(1,3,5), qstep=50e6,
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, 
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
dev.off()


##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PLOT HYGIDA
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


paf <- bcop_bhyg_paf
head(paf)
main<- "Bradysia hygida to Bradysia coprophila"; 
ylabels <- c("Bcop", "Bhyg")
ttigcol = rgb(0.267004, 0.004874, 0.329415, 1);
qtigcol = rgb(0.28229, 0.145912, 0.46151, 1.0)
chrnames <- c("X","II","III","IV", "chrA", "chrB", "chrC", "chrX")
tnames = unique(paf[order(paf$tlen, decreasing = TRUE),6:7])$target ; head(tnames) ; length(tnames)
tassocnames = tnames[!(tnames %in% c("X","II","III","IV"))] ; head(tassocnames) ; length(tassocnames)
torder = c("X","II","III","IV", tassocnames) ; torder
qnames = unique(paf[order(paf$qlen, decreasing = TRUE),1:2])$query ; head(qnames) ; length(qnames)
qassocnames = qnames[!(qnames %in% c("chrX","chrA", "chrB", "chrC"))] ; head(qassocnames) ; length(qassocnames)
qorder = c("chrX","chrA", "chrB", "chrC", qassocnames) ; qorder
qchr <- paf$query %in% chrnames ; sum(qchr)
tchr <- paf$target %in% chrnames ; sum(tchr)
chrpaf <- paf[qchr & tchr,] ; dim(chrpaf)
bhygchrpaf <- chrpaf
head(bhygchrpaf)

pdf("Rfigures/bcop_bhyg_dot.ALTstepsize.pdf", width = 7, height = 7)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=bhygchrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", qstep=90e6, #grid.subset.q = c(1,3,5),
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, 
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
#axis(side=2,at=seq(0,600e6,100e6), labels = seq(0,600,100), las=1, font=2)
dev.off()



##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##
## PLOT SELF SELF
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


paf <- bcop_bcop_paf
head(paf)
main<- "Bradysia coprophila to Bradysia coprophila"; 
ylabels <- c("Bcop", "Bcop")
ttigcol = rgb(0.267004, 0.004874, 0.329415, 1);
qtigcol = rgb(0.267004, 0.004874, 0.329415, 1)
chrnames <- c("X","II","III","IV")
tnames = unique(paf[order(paf$tlen, decreasing = TRUE),6:7])$target 
tassocnames = tnames[!(tnames %in% c("X","II","III","IV"))] 
torder = c("X","II","III","IV", tassocnames) ; torder
qnames = tnames 
qassocnames = tassocnames
qorder = torder
qchr <- paf$query %in% chrnames ; sum(qchr)
tchr <- paf$target %in% chrnames ; sum(tchr)
chrpaf <- paf[qchr & tchr,] ; dim(chrpaf)
bhygchrpaf <- chrpaf


pdf("Rfigures/bcop_bcop_dot.pdf", width = 7, height = 7)
par(mar=c(6.1,6.1,3.1,2.1))
assembalign.dotplot(paf=bhygchrpaf, targetOrder = torder, queryOrder = qorder, 
                    grid.lines = TRUE, grid.lwd = 1, grid.col = "black", 
                    qspecificticks = TRUE, pos.goncol=goncol, neg.goncol = ngoncol, bgoncol=bgoncol, 
                    qtigcol = qtigcol, ttigcol = ttigcol, querygaps = qgaps, targetgaps = tgaps, 
                    gapcol=gapcol, gapbcol=gapcol, plotqueryticks = TRUE, xticknorm.target = 1e6, 
                    xticknorm.query = 1e6, xlab = "Pos (Mb)", querytiglabels = TRUE, 
                    targettiglabels = TRUE, font=2, font.lab=2, font.axis=2, ylabels = ylabels, 
                    xtext.line = 2.5, xtext.cex = 1, bty="n", tspecificticks = TRUE,
                    segwd=circsize)
dev.off()




