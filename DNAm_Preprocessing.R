




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


## set working directory to folder with all idats (96 files=48 samples)
#setwd("F:/Cluster/Collaborations/MitoLab (Martin)/Project_PIC_12958_B01_CUS_MethylEPIC.2017-12-11/idat/idats_all")
#baseDir<-"F:/Cluster/Collaborations/MitoLab (Martin)/Project_PIC_12958_B01_CUS_MethylEPIC.2017-12-11/idat/idats_all"
setwd("/Users/gabrielsturm/Google Drive (gabriel.sturm@nyspi.columbia.edu)/MitoLab - General/DATA/Data Core Labs/DNA Methylation/UCNG_2019-9189_Picard_Meth_Epic_Lifespan_Study/idats")
dir()
baseDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/DATA/Data Core Labs/DNA Methylation/UCNG_2019-9189_Picard_Meth_Epic_Lifespan_Study/idats/all_idats"

files <-list.files(path = baseDir)
length(files)
files

baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2) # Mac
dir()
basenames <- read.csv("basenames.csv")

### Generate Folder with all idat files
# setwd(baseDir)
# file_list <- list.files(path = baseDir, pattern = "*.idat", recursive = TRUE, full.names = TRUE)
# dir.create("all_idats")
# all_idats_dir <- "/Users/gabrielsturm/Google Drive (gabriel.sturm@nyspi.columbia.edu)/MitoLab - General/DATA/Data Core Labs/DNA Methylation/UCNG_2019-9189_Picard_Meth_Epic_Lifespan_Study/idats/all_idats/"
# for(i in 1:length(file_list)) {
#   file.copy(file_list[i], all_idats_dir)
# }

### Setup Targets File
# setwd("/Users/gabrielsturm/Google Drive (gabriel.sturm@nyspi.columbia.edu)/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 3- Epigenetic Cell Aging/DM_R_Analysis/")
# load("betas.rcp.FunNorm.rdata")
# write.csv(targets, "targets.csv")

## Read.rgSet

## Merge with experimental info
targets <- read.metharray.sheet(baseDir)
length(unique(targets$Basename))  ## 512 unique arrays

## Read Methylation data from IDAT files
Raw.RGset<-read.metharray.exp(targets = targets)

## Match with Targets
length(intersect(targets$Sample_Name, Raw.RGset$Sample_Name)) ## 512

identical(targets$Sample_Name,Raw.RGset$Sample_Name) # TRUE
all(targets$Sample_Name==Raw.RGset$Sample_Name) # TRUE

setwd(baseDir2)
save(Raw.RGset, targets, file = "")

load("RGset.Raw.RData")

# Get Mset
MSet <- preprocessRaw(Raw.RGset)
methylated <- getMeth(MSet)
methylated[1:5,1:5]
unmethylated <- getUnmeth(MSet)
unmethylated[1:5,1:5]

setwd(baseDir2)
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
setwd(baseDir2)
save(dataMethylAid,targets,file="dataMethylAid.RAW.RData")
visualize(dataMethylAid)




## Sex-prediction with genomic ranges for X-Y chromosomes
setwd(baseDir2)
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
baseDir3 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp"
setwd(baseDir3)
dir()
LS_Data <- read.csv("Lifespan_Study_Data.csv")

Bad_Sample_info <- LS_Data[LS_Data$basename %in% QC_Bad_Samples$x,]
nrow(Bad_Sample_info) == nrow(QC_Bad_Samples)
Bad_Sample_info$Unique_Variable_Name

Bad_Sample_info_2 <- LS_Data[LS_Data$basename %in% other_bad_samples,]
nrow(Bad_Sample_info_2) == length(other_bad_samples)
Bad_Sample_info_2$Unique_Variable_Name

## Generate MethylSet from Genomic Ratio Set
MSet.raw <- preprocessRaw(Raw.RGset)
save(MSet.raw,file="MSet_raw.RData")


setwd(baseDir)
targets <- read.csv("targets_2.csv", row.names = 1)

## Same dimensions?
dim(targets);dim(MSet.raw)
targets <- targets[match(colnames(MSet.raw), rownames(targets)),]
all(rownames(targets)==colnames(MSet.raw)) ## TRUE
identical(rownames(targets),colnames(MSet.raw)) ## TRUE

setwd(baseDir2)
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
##
##
##   Detecion P-values
##
##
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
##                 After Fun Norm
##
##
##
##
##
##
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
##
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
##
##
##
##
##
##

library(sva)                                           # ComBat Adjustment

setwd(baseDir2)
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


setwd(baseDir2)
save(betas, targets, file = "Combat_Betas.RData")

# make GEO processed matrix
setwd(baseDir2)
load("DetectionPvalues.RData")
detP[1:5,1:5]
load("Combat_Betas.RData")
dim(betas)
detP <- detP[rownames(detP) %in% rownames(betas),]
nrow(detP) == nrow(betas)
detP <- detP[match(rownames(betas),rownames(detP)),]
identical(rownames(betas),rownames(detP))

# Exclude MiSBIE Samples
targets_2 <- subset(targets, targets$Experiment != 5)
nrow(targets_2) # 496
betas <- betas[,colnames(betas) %in% rownames(targets_2)]
ncol(betas) # 496

# remove excluded samples
bad_samples <- unique(c(QC_Bad_Samples$x, other_bad_samples,"203784950019_R01C01"))
length(bad_samples)
betas <- betas[,!colnames(betas) %in% bad_samples]
ncol(betas) # 480

# match samples in detP matrix
detP <- detP[,colnames(detP) %in% colnames(betas)]
detP <- detP[,match(colnames(betas),colnames(detP))]
identical(colnames(betas),colnames(detP))

# combine matrices
rows.combined <- nrow(betas) 
cols.combined <- ncol(betas) + ncol(detP)
GEO_processed_df <- matrix(NA, nrow=rows.combined, ncol=cols.combined)
GEO_processed_df[, seq(1, cols.combined, 2)] <- betas
GEO_processed_df[, seq(2, cols.combined, 2)] <- detP
rownames(GEO_processed_df) <- rownames(betas)
colnames(GEO_processed_df) <- paste0(rep(colnames(betas),each=2),"",rep(c(" beta"," Detection Pval"), ncol(betas)))
GEO_processed_df[1:5,1:10]
dim(GEO_processed_df)

setwd(baseDir2)
write.csv(GEO_processed_df,"GEO_DNAm_processed_matrix.csv")


# make signal intensity matrix for GEO submission
setwd(baseDir2)
load("DetectionPvalues.RData")
detP[1:5,1:5]
load("Methylated_Signal.RData")
load("Unmethylated_Signal.RData")
nrow(detP) == nrow(unmethylated)
identical(rownames(detP), rownames(unmethylated))
nrow(detP) == nrow(methylated)
identical(rownames(detP), rownames(methylated))

# Exclude MiSBIE Samples
targets_2 <- subset(targets, targets$Experiment != 5)
nrow(targets_2) # 496
methylated <- methylated[,colnames(methylated) %in% rownames(targets_2)]
ncol(methylated) # 496

# remove excluded samples
bad_samples <- unique(c(QC_Bad_Samples$x, other_bad_samples,"203784950019_R01C01"))
length(bad_samples)
methylated <- methylated[,!colnames(methylated) %in% bad_samples]
ncol(methylated) # 480

# match samples in detP matrix
detP <- detP[,colnames(detP) %in% colnames(methylated)]
detP <- detP[,match(colnames(methylated),colnames(detP))]
identical(colnames(methylated),colnames(detP))

# match samples in unmethylated matrix
unmethylated <- unmethylated[,colnames(unmethylated) %in% colnames(methylated)]
unmethylated <- unmethylated[,match(colnames(methylated),colnames(unmethylated))]
identical(colnames(methylated),colnames(unmethylated))

# combine matrices
rows.combined <- nrow(unmethylated) 
cols.combined <- ncol(unmethylated) + ncol(methylated) + ncol(detP)
GEO_signal_df <- matrix(NA, nrow=rows.combined, ncol=cols.combined)
GEO_signal_df[, seq(1, cols.combined, 3)] <- unmethylated
GEO_signal_df[, seq(2, cols.combined, 3)] <- methylated
GEO_signal_df[, seq(3, cols.combined, 3)] <- detP
rownames(GEO_signal_df) <- rownames(unmethylated)
colnames(GEO_signal_df) <- paste0(rep(colnames(methylated),each=3)," ",rep(c("Unmethylated Signal","Methylated Signal","Detection Pval"), ncol(methylated)))
GEO_signal_df[1:5,1:10]
dim(GEO_signal_df)
setwd(baseDir2)
write.csv(GEO_signal_df,"GEO_DNAm_signal_intensities_matrix.csv")


# metadata for GEO sumbission
setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp")
dir()
GEO_metadata <- read.csv("downloadable_data/Cellular_Lifespan_study_DNAm.csv")
GEO_metadata <- GEO_metadata[!is.na(GEO_metadata$DNAmethylation_sampleID),]
GEO_metadata <- GEO_metadata[GEO_metadata$DNAmethylation_sampleID != "",]
GEO_metadata <- GEO_metadata[order(GEO_metadata$DNAmethylation_sampleID),]
rownames(GEO_metadata) <- GEO_metadata$DNAmethylation_sampleID
nrow(GEO_metadata)

# remove excluded samples
bad_samples <- unique(c(QC_Bad_Samples$x, other_bad_samples,"203784950019_R01C01"))
length(bad_samples)
GEO_metadata <- GEO_metadata[!GEO_metadata$DNAmethylation_sampleID %in% bad_samples,]
nrow(GEO_metadata) # 479

# keep rows in sensical order
GEO_metadata <- GEO_metadata[order(GEO_metadata$Sample),]

# # find two missing samples
targets_3 <- targets_2[rownames(targets_2) %in% colnames(methylated),]
missing_metadata_samples <- targets_3[!rownames(targets_3) %in% GEO_metadata$DNAmethylation_sampleID,]
nrow(missing_metadata_samples) # 0

# add info from targets file
targets_3 <- targets_3[match(rownames(GEO_metadata),rownames(targets_3)),]
identical(rownames(GEO_metadata),rownames(targets_3))
GEO_metadata <- data.frame(GEO_metadata,targets_3)

setwd(baseDir2)
write.csv(GEO_metadata,"GEO_DNAm_metadata.csv")









## Run PCA on entire dataset
setwd(baseDir2)
dir()
load("betas.rcp.FunNorm.RData")

# Exclude MiSBIE Samples
targets_2 <- subset(targets, targets$Experiment != 5)
nrow(targets_2) # 496
data <- beta.rcp[,colnames(beta.rcp) %in% rownames(targets_2)]
ncol(data) # 496

# Exclude HEK Cells
targets_3 <- subset(targets_2, targets_2$person != "HEK293")
nrow(targets_3) # 484
data <- data[,colnames(data) %in% rownames(targets_3)]
ncol(data) # 484

pca = prcomp(t(data), center = TRUE, scale = TRUE) # 10 min
print(pca)
sum <- summary(pca)
sum$importance[2] # Proportion of Variance Explained for PC1
sum$importance[5] # Proportion of Variance Explained for PC2
sum$importance[8] # Proportion of Variance Explained for PC3

plot(pca, type = "l", ylim=c(0, 150000))
PoV <- pca$sdev^2/sum(pca$sdev^2)*100
plot(PoV, type = "o", ylab = "Percent of Variance", xlab = "PC", ylim = c(0,20),xlim= c(1,10), xaxp  = c(0,10,10))

pca$x[,1:3]

library(pca3d)
pca2d(pca, group = colnames(data), show.labels = FALSE)
axis(1, xaxp = c(-900,900,20))
axis(2, yaxp = c(-800,500,20))
pca3d(pca, group = colnames(data), show.labels = FALSE)

#export PC's 1-3
write.csv(pca$x[,1:5], "PCA_Data_HEK_removed.csv")

targets_3 <- targets_3[match(rownames(pca$x), rownames(targets_3)),]
targets_3$group <- paste(targets_3$person,targets_3$status, sep = "_")
pca_data <- cbind(targets_3, pca$x[,1:3])

shapes <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)
ggplot(pca_data, aes(x=PC1, y =PC2, color = person), label = age) +
  geom_point(size = 4 , stroke = 1, aes(shape = status)) +
  #geom_text(aes(label=Percent_Change),hjust=-0.2, vjust=1.2, size = 3) +
  #scale_color_manual(values = colors) +
  #scale_fill_manual(values = colors_2) +
  geom_segment(aes(
    xend=c(tail(PC1, n=-1), NA), 
    yend=c(tail(PC2, n=-1), NA)), alpha = 0.1) + 
  scale_shape_manual(values=shapes) +
  geom_text(aes(label=paste(round(age, digits = 0), "days", sep = " ")),hjust=-.4, vjust=-.4, size = 3, alpha = 0.7) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  #scale_y_continuous(name = Y_title) +
  #scale_x_continuous(name =  X_title) +
  theme_classic() +
  #annotate("text", label = annotation, x = 2, y = 0, size = 3, hjust = 0) +
  theme(text = element_text(size = 14),
        #legend.position="none",
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        plot.margin = margin(t = 6, r = 6, b = 6, l = 6, unit = "pt"),
        plot.title = element_text(size = 32, hjust = 0.05, vjust = -.1))


# Save beta.rcp with new targets file
setwd(baseDir2) # Mac
targets <- read_excel("Lifespan_DNAm_Samples.xlsx") %>% janitor::clean_names()
save(beta.rcp,targets,file="betas.rcp.FunNorm.RData", compress = T)
