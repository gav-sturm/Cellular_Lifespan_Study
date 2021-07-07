########
# Title: DNA Methylation Preprocessing Script
# Author: Gabriel Sturm (adapted from Andres Cardenas)
# Date: 2020-01-01
# see MinFi pipeline: https://www.bioconductor.org/packages/release/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
#

## Load packages needed
R.version$version.string  

## BioC packages needed
methpackagesBioC <- c("IlluminaHumanMethylation450kanno.ilmn12.hg19",
                      "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
                      "minfi",
                      "bumphunter",
                      "IlluminaHumanMethylationEPICmanifest",
                      "sva")

## Install from BioC
toinstallBioC <- setdiff(methpackagesBioC, installed.packages()[,1])
if(length(toinstallBioC >= 1)) {
  source("https://bioconductor.org/biocLite.R")
  biocLite(toinstallBioC, suppressUpdates = T)
  cat("finished installing new packages from BioConductor\n")
} else cat("packages we need from BioConductor are already installed\n")


## Packcages needed
library(minfi)
library(IlluminaHumanMethylationEPICmanifest)


## set directory to folder with all idats (2 files per sample)
rawFileDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/DATA/Data Core Labs/DNA Methylation/UCNG_2019-9189_Picard_Meth_Epic_Lifespan_Study/idats/all_idats"
setwd(rawFileDir)
files <-list.files(path = rawFileDir)
length(files)
files

## set working directory to place processed dataframes
workingDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(workingDir) # Mac
dir()
basenames <- read.csv("basenames.csv")


## Read.rgSet

## Merge with experimental info
targets <- read.metharray.sheet(rawFileDir)
length(unique(targets$Basename))  ## 512 unique arrays

## Read Methylation data from IDAT files
Raw.RGset<-read.metharray.exp(targets = targets)

## Match with Targets
length(intersect(targets$Sample_Name, Raw.RGset$Sample_Name)) ## 512

identical(targets$Sample_Name,Raw.RGset$Sample_Name) # TRUE
all(targets$Sample_Name==Raw.RGset$Sample_Name) # TRUE

setwd(workingDir)
save(Raw.RGset, targets, file = "")

load("RGset.Raw.RData")

# Get Mset
MSet <- preprocessRaw(Raw.RGset)
methylated <- getMeth(MSet)
methylated[1:5,1:5]
unmethylated <- getUnmeth(MSet)
unmethylated[1:5,1:5]

# save methylated and unmethylated signal matrices
setwd(workingDir)
save(methylated, file="Methylated_Signal.RData")
save(unmethylated, file="Unmethylated_Signal.RData")


## Examine using shiny
source("https://bioconductor.org/biocLite.R")
biocLite("shinyMethyl")
require(shinyMethyl)
myShinyMethylSet <- shinySummarize(Raw.RGset)

## Launch ShinyMethyl
runShinyMethyl(myShinyMethylSet)

## QC using methylAID
library(MethylAid)
library(BiocParallel)
require(minfi)
dataMethylAid <- summarize(targets)

## Save MethylAid objects
setwd(workingDir)
save(dataMethylAid,targets,file="dataMethylAid.RAW.RData")
visualize(dataMethylAid)


## Sex-prediction with genomic ranges for X-Y chromosomes
setwd(workingDir)
load("RGset.Raw.RData")
GRset <- mapToGenome(Raw.RGset)

## Plot sex-prediction for the entire dataset
tiff("Sex_Predicted.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

identical(Raw.RGset$Sample_Name,targets$Sample_Name)
all(Raw.RGset$Sample_Name==targets$Sample_Name)

plotSex(getSex(GRset, cutoff = -2),id=targets$Sample_Name)

dev.off()

sex_prediction <- getSex(GRset)

# save sex predictions
write.csv(sex_prediction, "DNAm_sample_sex_prediction.csv")



## QCplot
qc<-getQC(GRset)
addQC(GRset, qc)
plotQC(qc, badSampleCutoff = 11)

## Make QC function
all(rownames(qc)==targets$Basename)
my.plotQC<-function (qc, badSampleCutoff = 0.97) 
{
  meds <- (qc$uMed/qc$mMed) # cutoff should be 10.5
  whichBad <- which((meds < badSampleCutoff))
  Color = c(ifelse(targets$Experiment==5,"deeppink","deepskyblue"), ifelse(targets$Experiment==4,"deeppink","green"))
  Symbol = ifelse(targets$Experiment==5,22,1)
  plot(qc$mMed, qc$uMed, xlim = c(10, 14), ylim = c(10, 14), pch=Symbol,
       xaxt = "n", yaxt = "n", xlab = "Meth median intensity (log2)", 
       ylab = "Unmeth median intensity (log2)", col=Color)
  axis(side = 1, at = c(10, 12, 14))
  axis(side = 2, at = c(10, 12, 14))
  abline(0.3, 0.95, lty = 2)
  if (length(whichBad) > 0) 
    text(qc$mMed[whichBad]+0.05, qc$uMed[whichBad]-0.05, labels = rownames(qc)[whichBad], 
         col = "red",pch = 4, cex = 0.4)
  legend("topleft", legend = c("bad, with sample index"), 
         pch = 16, col = "red", bty = "n")
  invisible(NULL)
  legend("topright", legend = unique(targets$Experiment), 
         pch =c(1,22), col = c("deepskyblue","deeppink","blue","green", "pink"), bty = "n")
  invisible(NULL)
}
## Densities by sample Plate
tiff("QC_plot.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

my.plotQC(qc)

dev.off()

# identify the bad samples below QC threshold
badSampleCutoff = 0.97
meds <- (qc$uMed/qc$mMed)
whichBad <- which((meds < badSampleCutoff))
QC_Bad_Samples <- rownames(qc)[whichBad]
write.csv(QC_Bad_Samples, "QC_Bad_Samples.csv")


### Find which are the bad samples
QC_Bad_Samples <- read.csv("QC_Bad_Samples.csv")
other_bad_samples <- c("203833030104_R02C01",
                       "203833030104_R03C01",
                       "203833030104_R04C01",
                       "203833030104_R08C01")
intersect(QC_Bad_Samples, other_bad_samples)

# load in full lifespan study datasheet
lifespanDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp"
setwd(lifespanDir)
LS_Data <- read.csv("Lifespan_Study_Data.csv")
setwd(workingDir)

Bad_Sample_info <- LS_Data[LS_Data$basename %in% QC_Bad_Samples$x,]
nrow(Bad_Sample_info) == nrow(QC_Bad_Samples)
Bad_Sample_info$Unique_Variable_Name

Bad_Sample_info_2 <- LS_Data[LS_Data$basename %in% other_bad_samples,]
nrow(Bad_Sample_info_2) == length(other_bad_samples)
Bad_Sample_info_2$Unique_Variable_Name

## Generate MethylSet from Genomic Ratio Set
MSet.raw <- preprocessRaw(Raw.RGset)
save(MSet.raw,file="MSet_raw.RData")


setwd(rawFileDir)
targets <- read.csv("targets_2.csv", row.names = 1)

## Same dimensions?
dim(targets);dim(MSet.raw)
targets <- targets[match(colnames(MSet.raw), rownames(targets)),]
all(rownames(targets)==colnames(MSet.raw)) ## TRUE
identical(rownames(targets),colnames(MSet.raw)) ## TRUE

setwd(workingDir)
## Density plot by sample plate
require(scales)
tiff("Sample_by_plate.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(MSet.raw, sampGroups=targets$Sample_Plate, pal =alpha(rainbow(length(unique(targets$Sample_Plate))),0.75)) #

dev.off()


## Density plot by Experiment
require(scales)
tiff("Sample_by_Experiment.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(MSet.raw, sampGroups=targets$Experiment, pal =alpha(c("red","blue","green","purple", "yellow"),0.75)) #

dev.off()

#install.packages("viridis")  # Install
library("viridis") 

## Density plot by DNeasy Experiment
require(scales)
tiff("Sample_by_DNeasy.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(MSet.raw, sampGroups=targets$Dneasy_Batch, pal =alpha(rainbow(length(unique(targets$Dneasy_Batch))),0.5)) #

dev.off()



#############################
## SNPs Check
##
##
## Get SNPs
SNPs_probes<-getSnpBeta(Raw.RGset)

## SNP as matrix
SNPs<-as.data.frame(SNPs_probes)

## Histograms of all of the SNPs

hist(as.matrix(SNPs), main="Histogram of 59 SNPs in 600 samples",col="grey88",
     xlab="Beta value from SNPs",ylab="Density",prob=TRUE,ylim = c(0,3))
lines(density(as.matrix(SNPs)),col="deepskyblue",lwd=2)

## Make sure tartgets align with Sentris
targets <- targets[match(colnames(SNPs), rownames(targets)),]
all(colnames(SNPs)==rownames(targets)) # TRUE
identical(colnames(SNPs),rownames(targets)) # TRUE


SNP.indicator <- as.data.frame(sapply(SNPs,function(x) 
  ifelse(x>=0.20 & x<=0.30,1,
         ifelse(x>=0.60 & x<=0.80,1,0))))


## Look at failed samples
identical(rownames(targets),colnames(SNP.indicator)) # TRUE

## Total.SNPs
total.SNPs<-colSums(SNP.indicator)
length(names(total.SNPs)) # 512

## Count of Mix-ups
require(Cairo)
require(RColorBrewer)
colors <- brewer.pal(10, "BuPu")


tiff("Sample_SNP.tiff", width = 24, height = 8, 
     units = 'in', res = 600, compression = "lzw")


## Graph detection P-values
barplot(total.SNPs, col=colors, las=3, 
        cex.names=0.45, ylab="",ylim=c(0,15))
abline(h=6,col="red")
title(ylab="#-SNPs with unclear beta range", line=5, cex.lab=1.2)

dev.off()


## Annotation of Infinium type for each probe (I vs II)
typeI <-   minfi::getProbeInfo(MSet.raw,type="I")$Name
typeII <-  minfi::getProbeInfo(MSet.raw,type="II")$Name
onetwo <- rep(1, nrow(MSet.raw))
onetwo[rownames(MSet.raw) %in% typeII] <- 2
# almost 84% of our probes are type II
table(onetwo)



## Density before FunNorm
require(scales)
tiff("Probe_type_All.tiff", width = 12, height = 8, 
     units = 'in', res = 400, compression = "lzw")


densityPlot(MSet.raw[rownames(getAnnotation(MSet.raw)) %in% typeI,],pal = alpha("red",0.30),main='Beta density')
densityPlot(MSet.raw[rownames(getAnnotation(MSet.raw)) %in% typeII,],add = F, pal = alpha("blue",0.30))
legend("topright", c("Infinium I","Infinium II"), 
       lty=c(1,1), title="Infinium type", 
       bty='n',col=c("red","blue"))

dev.off()



## Get QC from minfi
out<-minfiQC(MSet.raw) # 5 min

# 5 min
save(out,file="Minfi_QCFile_Raw.RData")     

# 5 min
qcReport(Raw.RGset, sampNames=targets$Sample_Name, sampGroups=targets$Sample_Plate,
         pdf="qcReport_LS_DNAm_2.pdf")



##################################################
##
##
##   Detecion P-values
##
##
##

detP <- detectionP(Raw.RGset)
head(detP)

save(detP,file="DetectionPvalues.RData")

## How many samples with P<0.001
keep <- colMeans(detP) < 0.001
table(keep)


## Check aligment of Pvals with targets
identical(colnames(detP),rownames(targets)) # TRUE
all(colnames(detP)==rownames(targets)) # TRUE
length(rownames(targets))==length(colnames(detP)) # TRUE

#colnames(detP)<-rownames(targets)



## Detection P-values
require(scales)
require(RColorBrewer)
tiff("Detection_Pvalues.tiff", width = 12, height = 8, 
     units = 'in', res = 400, compression = "lzw")
par(mfrow=c(1,1),mai = c(0.7, 1, 0.1, 0.1))


pal <- brewer.pal(8,"Dark2")

## Graph detection P-values
barplot(colMeans(detP), col=pal, las=2, 
        cex.names=0.7, ylim=c(0,0.02), ylab="P Value")
abline(h=0.001,col="red")
title(ylab="Mean detection p-values", line=5, cex.lab=1.2)



dev.off()

Pvalue_bad <- which(colMeans(detP) > 0.01)
write.csv(Pvalue_bad, "Detection_P_Values_Bad_Samples.csv")


################################################################
## 
##
##  Fun Norm
##
##
##
##



## Normalize Samples
Mset.norm<-preprocessFunnorm(Raw.RGset, nPCs=2, bgCorr = TRUE,dyeCorr = TRUE, verbose = TRUE, ratioConvert = F) # 20 min


identical(rownames(targets),colnames(Mset.norm)) # TRUE
all(rownames(targets)==colnames(Mset.norm)) # TRUE
class(Mset.norm)

## Density After FunNorm
require(scales)
tiff("Sample_by_PlateFunNorm.tiff", width = 12, height = 8, 
     units = 'in', res = 600, compression = "lzw")

densityPlot(getBeta(Mset.norm), sampGroups=targets$Sample_Plate, pal =alpha(c("red","blue"),0.75)) #

dev.off()

## Save FunNorm mSet
save(Mset.norm,targets,file = "Mset.FunNorm.RData")


##########################
## RCP
##
##

## Aligment
all(rownames(targets)==colnames(Mset.norm)) # TRUE
## Identical
identical(rownames(targets),colnames(Mset.norm)) # TRUE


## densityPlot does not take Mset.Quantile (Genomic-RatioSet)
require(ENmix)
require(minfi)

Mset.norm <- updateObject(Mset.norm)


## densityPlot does not take Mset.Quantile (Genomic-RatioSet)
betas.FunNorm <- minfi::getBeta(Mset.norm, type = "Illumina") # 1min
range(betas.FunNorm) # 0 to 0.99
median(betas.FunNorm) # 0.58


## Get values from FunNorm Data
raw.M <- logit2(betas.FunNorm) # 30 sec
range(raw.M) # 0 to 8.2
## Check aligment
all(rownames(targets)==colnames(raw.M)) # TRUE
identical(rownames(targets),colnames(raw.M)) # TRUE


## Other objects needed
dist=25
quantile.grid=seq(0.001,0.999,by=0.001)
qcscore=NULL
nbthre=3
detPthre=0.000001

# find therby pairs of type I probes and type II probes
GRset <- mapToGenome(Mset.norm)
annotation<-getAnnotation(GRset) # 5 min
annotation=annotation[intersect(rownames(betas.FunNorm),rownames(annotation)),]
dim(annotation);dim(betas.FunNorm)

probe.II.Name = annotation$Name[annotation$Type=="II"]
annotation = annotation[order(annotation$chr,annotation$pos),]
anno1 = annotation[1:(nrow(annotation)-1),]
anno2 = annotation[2:nrow(annotation),]
flag = (abs(anno1$pos-anno2$pos)<dist & anno1$chr==anno2$chr & 
        anno1$Relation_to_Island==anno2$Relation_to_Island & anno1$Type !=
        anno2$Type)
anno1 = anno1[flag,]
anno2 = anno2[flag,]
probe.I = anno1$Name
probe.II = anno2$Name
probe.I[anno2$Type=="I"] = anno2$Name[anno2$Type=="I"]
probe.II[anno1$Type=="II"] = anno1$Name[anno1$Type=="II"]

raw.M.t = raw.M[c(probe.I,probe.II),]

#remove low quality data
if(is.null(qcscore)){}else if((sum(!(rownames(raw.M.t) %in% 
                                     rownames(qcscore$detP))) +
                               sum(!(colnames(raw.M.t) %in% colnames(qcscore$detP))))>0){
  stop("Wrong qcscore matrix, please check...\n")}else{
    temp <- qcscore$nbead<nbthre | qcscore$detP>detPthre
    temp=temp[rownames(raw.M.t),]
    temp=temp[,colnames(raw.M.t)]
    raw.M.t[temp]=NA
  }
# NUll

#linear regression
M.II <- raw.M.t[probe.II,]
M.I <- raw.M.t[probe.I,]


##quantile.grid=seq(0.001,0.999,by=0.001)
quantile.grid = seq(0.001,0.999,by=0.001)
qtl <- function(x) quantile(x, quantile.grid, na.rm=TRUE)
M.I = apply(M.I,2,qtl)
M.II = apply(M.II,2,qtl)

beta.est <- mat.or.vec(2,ncol(betas.FunNorm))

for (i in 1:ncol(betas.FunNorm)){
  index<-(M.II[,i]!=Inf & M.II[,i]!=-Inf & M.I[,i]!=Inf & M.I[,i]!=-Inf)
  X<-cbind(rep(1,sum(index)),M.II[index,i]); Y<-M.I[index,i]
  beta.est[,i]<-solve(t(X)%*%X)%*%t(X)%*%Y
}

M.II.all <- raw.M[probe.II.Name,]
M.II.new <- mat.or.vec(nrow(M.II.all),ncol(M.II.all))
for (i in 1:ncol(M.II.all)){
  M.II.new[,i]<-beta.est[1,i]+beta.est[2,i]*M.II.all[,i]
}
M.II.new[M.II.all==Inf]<-Inf; M.II.new[M.II.all==-Inf]<-(-Inf)

betas.FunNorm[probe.II.Name,] <- ilogit2(M.II.new) 
beta.rcp <- betas.FunNorm
range(beta.rcp) # 0 to 0.99
sum(is.na(beta.rcp)) # 0


## Aligment with pDat
all(colnames(beta.rcp)==rownames(targets)) # TRUE
identical(colnames(beta.rcp),rownames(targets)) # TRUE

save(beta.rcp,targets,file="betas.rcp.FunNorm.RData", compress = T)


typeI <-   minfi::getProbeInfo(Raw.RGset,type="I")$Name
typeII <-  minfi::getProbeInfo(Raw.RGset,type="II")$Name
onetwo <- rep(1, dim(Mset.norm)[1])
onetwo[rownames(Mset.norm) %in% typeII] <- 2
table(onetwo)

identical(colnames(Mset.norm),colnames(beta.rcp))
identical(rownames(Mset.norm),rownames(beta.rcp))

## Density by Infinium-Type
## Density After FunNorm
require(scales)
tiff("FunNorm_RCP.tiff", width = 12, height = 8, 
     units = 'in', res = 400, compression = "lzw")



densityPlot(beta.rcp[rownames(beta.rcp) %in% typeI,],pal = alpha("red",0.65),main='Beta density')
densityPlot(beta.rcp[rownames(beta.rcp) %in% typeII,],add = F, pal = alpha("blue",0.65))
legend("topright", legend = c("Infinium I", "Infinium II"), 
       lty =c(1,1), col = c("red","blue"), bty = "n")

dev.off()



###############################
## Combat Adjustment
##
##

library(sva)                                           # ComBat Adjustment

setwd(workingDir)
load("betas.rcp.FunNorm.RData")

setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/")
targets <- read.csv("targets_2.csv", row.names = 1)
table(targets$Sample_Plate) # 96 by 96 by 96 by 96 by 96 by 32
#' Adjust for Batch=Plate
length(unique(targets$Sample_Plate))                   # 6-plates
sum(is.na(beta.rcp)) # 0 missing betas
range(beta.rcp) # 0 to 0.999

mvals <-log2(beta.rcp)-log2(1-beta.rcp) # transform to M-values
range(mvals) # 0 to 9
sum(is.na(mvals)) # 0

#' ComBat Adjust Data
## RowVariance can't be 0
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE,
                    SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}


# Calculate variance
vars = as.matrix(rowVars(mvals))

# Replace all probes with no variance with NA
# and remove them from the FunNorm set
vars[vars == 0] = NA
vars = na.omit(vars)

## How many NAs
sum(is.na(vars)) # 0
mval.filter <- mvals[rownames(mvals) %in% rownames(vars),]
dim(mval.filter);dim(mvals) # missing ~200 sites

## Make sure there is aligment before running
targets <- targets[match(colnames(mvals.ComBat), rownames(targets)),]
identical(rownames(targets), colnames(mval.filter)) # TRUE
all(rownames(targets)==colnames(mval.filter)) # TRUE


# ComBat model
#table(as.numeric(as.character(targets$age)))
#targets$age<-as.numeric(as.character(targets$Chronological.Age..days.))
model <- model.matrix(~as.numeric(targets$age))

## This is the statistical model
model

# Batch
length(unique(targets$Sample_Plate)) # 6 plates
# Plates <- as.numeric(unique(targets$Sample_Plate))
# targets$Batch<- 
# table(targets$Batch,targets$Sample_Plate)

## ComBat
mvals.ComBat = ComBat(mval.filter,batch = targets$Sample_Plate, mod = model)

## check dimentions of cov
dim(model)
dim(mvals.ComBat)

## Check Aligment
identical(colnames(mvals.ComBat), rownames(targets)) # TRUE

length(unique(targets$Slide))                   # 64-slides

## Model Matrix Adjusted Mercury
Z = model.matrix(~as.numeric(targets$age)+    ## Main effect of Age 
                   as.factor(targets$Slide))   ## Let's adjust for slide

##Point of interest
vPOI = 2     ## THIS IS THE POINT OF INTEREST
table(Z[,2]) ## This is age

## Transform back to Beta values (i.e. undo logit transformation)
expit2 = function(x) 2^x/(1+2^x)
betas = expit2(mvals.ComBat)
range(betas) # 0 to 0.999
hist(betas) # grapg

## Check Aligment
identical(colnames(betas), rownames(targets)) # TRUE

# Save final Combat-adjusted dataframe
setwd(workingDir)
save(betas, targets, file = "Combat_Betas.RData")
