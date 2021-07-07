# rm(list = ls())
# gc()
# .rs.restartR()
library(foreach)
library(doParallel)
library(parallel)
library(lme4)
library(MuMIn)
library(ggplot2)
library(gridExtra)
library(grid)
library(DescTools)
library(rlist)
library(tidyverse)
library(qlcMatrix)
library("illuminaHumanv3.db") #BiocManager::install("illuminaHumanv3.db")
library(ggpmisc)
library(GEOquery)
library(robustbase)
library(miceadds)

baseDir<- c("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/Modeling_Chris_Kempes/")
setwd(baseDir)
dir()

######## PART 0 - Data Wrangling ########

# Lifespan Metadata
setwd(baseDir)
LS_Data <- read.csv('Cellular_Lifespan/Lifespan_Study_Data.csv')


##### Fibroblasts Controls r1 DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
dim(betas)
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

 # Save Data
setwd(baseDir)
control_fibroblasts_r1_betas <- filtered_fbr_betas
control_fibroblasts_r1_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r1_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_betas.RData")
write.csv(control_fibroblasts_r1_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_metadata.csv")
#######

##### Fibroblasts Controls r1 early-phase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 15 & filtered_fbr_metadata$Days_Grown < 60,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r1_earlyPhase_betas <- filtered_fbr_betas
control_fibroblasts_r1_earlyPhase_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r1_earlyPhase_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_betas.RData")
write.csv(control_fibroblasts_r1_earlyPhase_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_metadata.csv")
#######

##### Fibroblasts Controls r1 mid-phase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown >= 60 & filtered_fbr_metadata$Days_Grown < 110,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r1_midPhase_betas <- filtered_fbr_betas
control_fibroblasts_r1_midPhase_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r1_midPhase_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_betas.RData")
write.csv(control_fibroblasts_r1_midPhase_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_metadata.csv")
#######

##### Fibroblasts Controls r1 senescent-phase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown >= 110,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r1_senescentPhase_betas <- filtered_fbr_betas
control_fibroblasts_r1_senescentPhase_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r1_senescentPhase_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_betas.RData")
write.csv(control_fibroblasts_r1_senescentPhase_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_metadata.csv")
#######

##### Fibroblasts Controls r2 DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB11","hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(3),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r2_betas <- filtered_fbr_betas
control_fibroblasts_r2_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r2_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_betas.RData")
write.csv(control_fibroblasts_r2_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_metadata.csv")
#######

##### Fibroblasts Controls r3 DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(4),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r3_betas <- filtered_fbr_betas
control_fibroblasts_r3_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r3_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_betas.RData")
write.csv(control_fibroblasts_r3_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_metadata.csv")
#######


##### Fibroblasts Contact Inhibition DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Contact_Inhibition_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(3),]
filtered_fbr_metadata <-  filtered_fbr_metadata[!is.na(filtered_fbr_metadata$Time_Point),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Time_Point > 0,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Contact Inhibition DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
save(filtered_fbr_betas,file="Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_betas.RData")
write.csv(filtered_fbr_metadata, "Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_metadata.csv")

#########

##### Fibroblasts Contact Inhibition earlyPhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Contact_Inhibition_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(3),]
filtered_fbr_metadata <-  filtered_fbr_metadata[!is.na(filtered_fbr_metadata$Time_Point),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Time_Point > 20 & filtered_fbr_metadata$Time_Point < 100,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
#filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 10 & filtered_fbr_metadata$Days_Grown < 70,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Contact Inhibition DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
save(filtered_fbr_betas,file="Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_betas.RData")
write.csv(filtered_fbr_metadata, "Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_metadata.csv")

#########

##### Fibroblasts Contact Inhibition latePhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Contact_Inhibition_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(3),]
filtered_fbr_metadata <-  filtered_fbr_metadata[!is.na(filtered_fbr_metadata$Time_Point),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Time_Point > 80,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
#filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown >= 70,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Contact Inhibition latePhase DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
save(filtered_fbr_betas,file="Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_betas.RData")
write.csv(filtered_fbr_metadata, "Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_metadata.csv")

#########

##### Fibroblasts Hypoxia DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_3",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2,3),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
min(filtered_fbr_metadata$Days_Grown)
max(filtered_fbr_metadata$Days_Grown)

# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Hypoxia DNAm Histogram", xlab = "beta", breaks = 100)
hypoxia_betas <- filtered_fbr_betas
hypoxia_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(hypoxia_betas,file="Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_betas.RData")
write.csv(hypoxia_DNAm_metadata, "Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_metadata.csv")
#########

##### Fibroblasts DEX earlyPhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == c("DEX_21"),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 25 & filtered_fbr_metadata$Days_Grown < 95,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts MitoModulators DNAm Histogram", xlab = "beta", breaks = 100)
DEX_betas <- filtered_fbr_betas
DEX_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(DEX_betas,file="Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_betas.RData")
write.csv(DEX_DNAm_metadata, "Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_metadata.csv")
#########

##### Fibroblasts Controls r1-treatments DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == "Control_21",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(1,2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 25 & filtered_fbr_metadata$Days_Grown < 95,]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)

# Save Data
setwd(baseDir)
control_fibroblasts_r1_earlyPhase_betas <- filtered_fbr_betas
control_fibroblasts_r1_earlyPhase_metadata <- filtered_fbr_metadata
save(control_fibroblasts_r1_earlyPhase_betas,file="Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_betas.RData")
write.csv(control_fibroblasts_r1_earlyPhase_metadata,"Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_metadata.csv")
#######

##### Fibroblasts Modulators earlyPhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == c("Modulators_21"),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 25 & filtered_fbr_metadata$Days_Grown < 95,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts MitoModulators DNAm Histogram", xlab = "beta", breaks = 100)
modulators_betas <- filtered_fbr_betas
modulators_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(modulators_betas,file="Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_betas.RData")
write.csv(modulators_DNAm_metadata, "Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_metadata.csv")
#########

##### Fibroblasts Oligo earlyPhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == c("Oligomycin_21"),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 25 & filtered_fbr_metadata$Days_Grown < 95,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts Oligomycin DNAm Histogram", xlab = "beta", breaks = 100)
oligo_betas <- filtered_fbr_betas
oligo_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(oligo_betas,file="Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_betas.RData")
write.csv(oligo_DNAm_metadata, "Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_metadata.csv")
#########

##### Fibroblasts SURF1 earlyPhase DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB6","hFB7","hFB8")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == c("Control_21"),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(2),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Days_Grown > 25 & filtered_fbr_metadata$Days_Grown < 95,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts SURF1 DNAm Histogram", xlab = "beta", breaks = 100)
SURF1_betas <- filtered_fbr_betas
SURF1_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(SURF1_betas,file="Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_betas.RData")
write.csv(SURF1_DNAm_metadata, "Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_metadata.csv")
#########

##### Fibroblasts 2-Deoxyglucose DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Treatment == c("2-Deoxyglucose_21"),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$Study_Part %in% c(4),]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="Fibroblasts 2-Deoxyglucose_ DNAm Histogram", xlab = "beta", breaks = 100)
deoxyglucose_betas <- filtered_fbr_betas
deoxyglucose_DNAm_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(deoxyglucose_betas,file="Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_betas.RData")
write.csv(deoxyglucose_DNAm_metadata, "Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_metadata.csv")
#########

##### HEK293 DNAm ######
# Load in fibroblast DNAm data
baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("HEK293")
filtered_fbr_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_fbr_metadata <-  filtered_fbr_metadata[filtered_fbr_metadata$basename != "",]
filtered_fbr_metadata$basename
nrow(filtered_fbr_metadata)
# Retrieve TPM for control and Oligo samples
filtered_fbr_betas <- betas[,colnames(betas) %in% as.vector(filtered_fbr_metadata$basename)]
# Order betas in same order as Lifespan Data
filtered_fbr_betas <- filtered_fbr_betas[,match(filtered_fbr_metadata$basename, colnames(filtered_fbr_betas))]
fibroblast_timepoints <- filtered_fbr_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints / 365
filtered_fbr_metadata <- data.frame(filtered_fbr_metadata, Timepoints_years)
## Check values match
all(colnames(filtered_fbr_betas)==filtered_fbr_metadata$basename)
identical(colnames(filtered_fbr_betas),as.character(filtered_fbr_metadata$basename))

filtered_fbr_betas[1:5,1:5]
dim(filtered_fbr_betas)
max(filtered_fbr_betas)
min(filtered_fbr_betas)
hist(as.matrix(filtered_fbr_betas), main="HEK293 DNAm Histogram", xlab = "beta", breaks = 100)
HEK293_betas <- filtered_fbr_betas
HEK293_metadata <- filtered_fbr_metadata
# Save Data
setwd(baseDir)
save(HEK293_betas,file="Cellular_Lifespan/DNAm/HEK293_betas.RData")
write.csv(HEK293_metadata, "Cellular_Lifespan/DNAm/HEK293_metadata.csv")
#######

##### TwinsUK Skin DNAm #######
# twinsUK_skin_DNAm_betas <- read.table(file="TwinsUK/DNAm_Skin/GSE90124_betasProcessed_Non-normalisedValues_322samples_2016.txt", row.names=1)
# colnames(twinsUK_skin_DNAm_betas) <- twinsUK_skin_DNAm_betas[1,]
# twinsUK_skin_DNAm_betas <- twinsUK_skin_DNAm_betas[-1,]
# #twinsUK_skin_DNAm_betas <- data.matrix(twinsUK_skin_DNAm_betas)
# twinsUK_skin_DNAm_betas <- as.matrix(sapply(twinsUK_skin_DNAm_betas, as.numeric)) 
# hist(twinsUK_skin_DNAm_betas, main="twinsUK Skin DNAm Histogram", xlab = "beta", breaks = 100)
# save(twinsUK_skin_DNAm_betas,file="TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData")
# library(GEOquery)
# twinsUK_skin_DNAm_matrix <- getGEO("GSE90124")
# save(twinsUK_skin_DNAm_matrix,file="TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_GSE90124_matrix.RData")
# twinsUK_skin_DNAm_betas <- twinsUK_skin_DNAm_matrix$GSE90124_series_matrix.txt.gz@assayData$exprs
# hist(twinsUK_skin_DNAm_betas, main="twinsUK Skin DNAm Histogram", xlab = "beta", breaks = 100)
# dim(twinsUK_skin_DNAm_betas)
# save(twinsUK_skin_DNAm_betas,file="TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData")

#######

##### TwinsUK Fat DNAm #######
#load('TwinsUK/DNAm_Fat/TwinsUK_FAT_DNAm_betas.RData')
twinsUK_fat_DNAm_betas <- read.table(file='TwinsUK/DNAm_Fat/MuTHER_Fat_450K_norm_AE_030913.txt',sep= "\t", row.names=1,header=T)
twinsUK_fat_DNAm_betas <- twinsUK_fat_DNAm_betas[-1,]
#twinsUK_fat_DNAm_betas <- as.matrix(sapply(twinsUK_fat_DNAm_betas, as.numeric)) 
#hist(as.matrix(twinsUK_fat_DNAm_betas), main="twinsUK Fat DNAm Histogram", xlab = "beta", breaks = 100)
save(twinsUK_fat_DNAm_betas,file="TwinsUK/DNAm_Fat/twinsUK_Fat_DNAm_betas.RData")
twinsUK_fat_DNAm_betas[1:5,1:5]
dim(twinsUK_fat_DNAm_betas)
#########

###### MESA DNAm #####
library(GEOquery)
# Paper: https://onlinelibrary.wiley.com/doi/full/10.1111/acel.13229
# Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56046
MESA_DNAm_matrix <- getGEO("GSE56046")
MESA_DNAm_metadata <- MESA_DNAm_matrix$GSE56046_series_matrix.txt.gz@phenoData@data
write.csv(MESA_DNAm_metadata, "MESA/MESA_DNAm_metadata.csv")
save(MESA_DNAm_matrix,file="MESA/MESA_DNAm_GSE90124_matrix.RData")
#########

###### MESA RNAarray #####
MESA_RNAseq_matrix <- getGEO("GSE56047")
MESA_RNAseq_metadata <- MESA_RNAseq_matrix$`GSE56047-GPL10558_series_matrix.txt.gz`@phenoData@data
MESA_RNAseq_metadata_2 <- MESA_RNAseq_matrix$`GSE56047-GPL13534_series_matrix.txt.gz`@phenoData@data
MESA_RNAseq_data <- MESA_RNAseq_matrix$`GSE56047-GPL10558_series_matrix.txt.gz`@assayData$exprs
save(MESA_RNAseq_data,file="MESA/RNAarray/MESA_RNAseq_data.RData")

#########


##### Fibroblasts Controls r1 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1-linear RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown < 90,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1-senescent RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 90,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1 earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 15 & filtered_metadata$Days_Grown < 60,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1 midPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown >= 60 & filtered_metadata$Days_Grown < 110,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1 senescentPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown >=  110,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r1 treatments RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r1_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r1_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r1_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_genes.RData")
write.csv(control_fibroblasts_r1_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r2 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB11","hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r2_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r2_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r2_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_genes.RData")
write.csv(control_fibroblasts_r2_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r3 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB13")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r3_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r3_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r3_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_genes.RData")
write.csv(control_fibroblasts_r3_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r3 earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB13")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 10 & filtered_metadata$Days_Grown < 75,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r3_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r3_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r3_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase_RNAseq_genes.RData")
write.csv(control_fibroblasts_r3_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls r3 earlyPhase13 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB13")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 10 & filtered_metadata$Days_Grown < 75,]
filtered_metadata <- filtered_metadata[filtered_metadata$Passage %in% c(5,7,9,11,13),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

control_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
control_fibroblasts_RNAseq_genes[1:5,1:5]
dim(control_fibroblasts_RNAseq_genes)
max(control_fibroblasts_RNAseq_genes)
min(control_fibroblasts_RNAseq_genes)
hist(as.matrix(control_fibroblasts_RNAseq_genes), main="Fibroblasts Controls RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
control_fibroblasts_r3_RNAseq_genes <- control_fibroblasts_RNAseq_genes
control_fibroblasts_r3_RNAseq_metadata <- filtered_metadata
save(control_fibroblasts_r3_RNAseq_genes,file="Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase13_RNAseq_genes.RData")
write.csv(control_fibroblasts_r3_RNAseq_metadata,"Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase13_RNAseq_metadata.csv")
#########

##### Fibroblasts Contact Inhibition RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB11")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Contact_Inhibition_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$Time_Point),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Time_Point > 0,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

contact_inhibition_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
contact_inhibition_fibroblasts_RNAseq_genes[1:5,1:5]
dim(contact_inhibition_fibroblasts_RNAseq_genes)
max(contact_inhibition_fibroblasts_RNAseq_genes)
min(contact_inhibition_fibroblasts_RNAseq_genes)
hist(as.matrix(contact_inhibition_fibroblasts_RNAseq_genes), main="Fibroblasts Contact Inhibition RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
contact_inhibition_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(contact_inhibition_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/contact_inhibition_fibroblasts_RNAseq_genes.RData")
write.csv(contact_inhibition_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/contact_inhibition_RNAseq_metadata.csv")
#########

##### Fibroblasts Contact Inhibition earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB11")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Contact_Inhibition_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$Time_Point),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Time_Point > 20 & filtered_metadata$Time_Point < 100,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes[1:5,1:5]
dim(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes)
max(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes)
min(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes)
hist(as.matrix(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes), main="Fibroblasts Contact Inhibition RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
contact_inhibition_earlyPhase_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(contact_inhibition_earlyPhase_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Contact Inhibition latePhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB11")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Contact_Inhibition_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2,3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$Time_Point),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Time_Point > 80,]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

contact_inhibition_latePhase_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
contact_inhibition_latePhase_fibroblasts_RNAseq_genes[1:5,1:5]
dim(contact_inhibition_latePhase_fibroblasts_RNAseq_genes)
max(contact_inhibition_latePhase_fibroblasts_RNAseq_genes)
min(contact_inhibition_latePhase_fibroblasts_RNAseq_genes)
hist(as.matrix(contact_inhibition_latePhase_fibroblasts_RNAseq_genes), main="Fibroblasts Contact Inhibition RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
contact_inhibition_latePhase_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(contact_inhibition_latePhase_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_fibroblasts_RNAseq_genes.RData")
write.csv(contact_inhibition_latePhase_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_RNAseq_metadata.csv")
#########

##### Fibroblasts Hypoxia RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB11")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatments == "Control",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Percent_Oxygen == 3,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2,3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
min(filtered_metadata$Days_Grown)
max(filtered_metadata$Days_Grown)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

hypoxia_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
hypoxia_fibroblasts_RNAseq_genes[1:5,1:5]
dim(hypoxia_fibroblasts_RNAseq_genes)
max(hypoxia_fibroblasts_RNAseq_genes)
min(hypoxia_fibroblasts_RNAseq_genes)
hist(as.matrix(hypoxia_fibroblasts_RNAseq_genes), main="Fibroblasts Hypoxia RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
hypoxia_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(hypoxia_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_genes.RData")
write.csv(hypoxia_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_metadata.csv")
#########


##### Fibroblasts Hypoxia + Contact Inhibition RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB11")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatments == "Contact_Inhibition",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Percent_Oxygen == 3,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2,3),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

hypoxia_contact_inhibition_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
hypoxia_contact_inhibition_fibroblasts_RNAseq_genes[1:5,1:5]
dim(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes)
max(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes)
min(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes)
hist(as.matrix(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes), main="Fibroblasts hypoxia_contact_inhibition RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_genes.RData")
write.csv(hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls DEX RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "DEX_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

dex_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
dex_fibroblasts_RNAseq_genes[1:5,1:5]
dim(dex_fibroblasts_RNAseq_genes)
max(dex_fibroblasts_RNAseq_genes)
min(dex_fibroblasts_RNAseq_genes)
hist(as.matrix(dex_fibroblasts_RNAseq_genes), main="Fibroblasts dexs RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
dex_fibroblasts_RNAseq_genes <- dex_fibroblasts_RNAseq_genes
dex_fibroblasts_RNAseq_metadata <- filtered_metadata
save(dex_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_genes.RData")
write.csv(dex_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Controls DEX earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "DEX_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(1,2),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

dex_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
dex_fibroblasts_RNAseq_genes[1:5,1:5]
dim(dex_fibroblasts_RNAseq_genes)
max(dex_fibroblasts_RNAseq_genes)
min(dex_fibroblasts_RNAseq_genes)
hist(as.matrix(dex_fibroblasts_RNAseq_genes), main="Fibroblasts dexs RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

setwd(baseDir)
dex_fibroblasts_RNAseq_genes <- dex_fibroblasts_RNAseq_genes
dex_fibroblasts_RNAseq_metadata <- filtered_metadata
save(dex_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(dex_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts SURF1 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB6","hFB7","hFB8")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

surf1_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
surf1_fibroblasts_RNAseq_genes[1:5,1:5]
dim(surf1_fibroblasts_RNAseq_genes)
max(surf1_fibroblasts_RNAseq_genes)
min(surf1_fibroblasts_RNAseq_genes)
hist(as.matrix(surf1_fibroblasts_RNAseq_genes), main="Fibroblasts SURF1 RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
surf1_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(surf1_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_genes.RData")
write.csv(surf1_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts SURF1 earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB6","hFB7","hFB8")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Control_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

surf1_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
surf1_fibroblasts_RNAseq_genes[1:5,1:5]
dim(surf1_fibroblasts_RNAseq_genes)
max(surf1_fibroblasts_RNAseq_genes)
min(surf1_fibroblasts_RNAseq_genes)
hist(as.matrix(surf1_fibroblasts_RNAseq_genes), main="Fibroblasts SURF1 RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
surf1_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(surf1_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(surf1_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Oligomycin RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Oligomycin_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

oligo_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
oligo_fibroblasts_RNAseq_genes[1:5,1:5]
dim(oligo_fibroblasts_RNAseq_genes)
max(oligo_fibroblasts_RNAseq_genes)
min(oligo_fibroblasts_RNAseq_genes)
hist(as.matrix(oligo_fibroblasts_RNAseq_genes), main="Fibroblasts oligo RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
oligo_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(oligo_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_genes.RData")
write.csv(oligo_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Oligomycin earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Oligomycin_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

oligomycin_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
oligomycin_fibroblasts_RNAseq_genes[1:5,1:5]
dim(oligomycin_fibroblasts_RNAseq_genes)
max(oligomycin_fibroblasts_RNAseq_genes)
min(oligomycin_fibroblasts_RNAseq_genes)
hist(as.matrix(oligomycin_fibroblasts_RNAseq_genes), main="Fibroblasts modulators RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
oligomycin_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(oligomycin_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(oligomycin_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Modulators RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Modulators_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

modulators_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
modulators_fibroblasts_RNAseq_genes[1:5,1:5]
dim(modulators_fibroblasts_RNAseq_genes)
max(modulators_fibroblasts_RNAseq_genes)
min(modulators_fibroblasts_RNAseq_genes)
hist(as.matrix(modulators_fibroblasts_RNAseq_genes), main="Fibroblasts modulators RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
modulators_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(modulators_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_genes.RData")
write.csv(modulators_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Modulators earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Modulators_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(2),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

modulators_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
modulators_fibroblasts_RNAseq_genes[1:5,1:5]
dim(modulators_fibroblasts_RNAseq_genes)
max(modulators_fibroblasts_RNAseq_genes)
min(modulators_fibroblasts_RNAseq_genes)
hist(as.matrix(modulators_fibroblasts_RNAseq_genes), main="Fibroblasts modulators RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
modulators_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(modulators_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(modulators_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts 2-deoxyglucose RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "2-Deoxyglucose_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

deoxyglucose_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
dim(deoxyglucose_fibroblasts_RNAseq_genes)
max(deoxyglucose_fibroblasts_RNAseq_genes)
min(deoxyglucose_fibroblasts_RNAseq_genes)
hist(as.matrix(deoxyglucose_fibroblasts_RNAseq_genes), main="Fibroblasts 2-deoxyglucose RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
deoxyglucose_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(deoxyglucose_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_genes.RData")
write.csv(deoxyglucose_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts 2-deoxyglucose earlyPhase RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "2-Deoxyglucose_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 25 & filtered_metadata$Days_Grown < 95,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

deoxyglucose_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
dim(deoxyglucose_fibroblasts_RNAseq_genes)
max(deoxyglucose_fibroblasts_RNAseq_genes)
min(deoxyglucose_fibroblasts_RNAseq_genes)
hist(as.matrix(deoxyglucose_fibroblasts_RNAseq_genes), main="Fibroblasts 2-deoxyglucose RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
deoxyglucose_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(deoxyglucose_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_genes.RData")
write.csv(deoxyglucose_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts 2-deoxyglucose earlyPhase-hFB13 RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB13")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "2-Deoxyglucose_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[filtered_metadata$Days_Grown > 10 & filtered_metadata$Days_Grown < 75,]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

deoxyglucose_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
dim(deoxyglucose_fibroblasts_RNAseq_genes)
max(deoxyglucose_fibroblasts_RNAseq_genes)
min(deoxyglucose_fibroblasts_RNAseq_genes)
hist(as.matrix(deoxyglucose_fibroblasts_RNAseq_genes), main="Fibroblasts 2-deoxyglucose RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
deoxyglucose_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(deoxyglucose_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase13_fibroblasts_RNAseq_genes.RData")
write.csv(deoxyglucose_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase13_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts betahydroxybutyrate RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "betahydroxybutyrate_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

betahydroxybutyrate_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
betahydroxybutyrate_fibroblasts_RNAseq_genes[1:5,1:5]
dim(betahydroxybutyrate_fibroblasts_RNAseq_genes)
max(betahydroxybutyrate_fibroblasts_RNAseq_genes)
min(betahydroxybutyrate_fibroblasts_RNAseq_genes)
hist(as.matrix(betahydroxybutyrate_fibroblasts_RNAseq_genes), main="Fibroblasts betahydroxybutyrate RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
betahydroxybutyrate_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(betahydroxybutyrate_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_genes.RData")
write.csv(betahydroxybutyrate_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_metadata.csv")
#########

##### Fibroblasts Galactose RNAseq ######
# Load in fibroblast RNAseq data
genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
fibroblast_rnaseq_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# Rename Samples
old_sample_names <- colnames(fibroblast_rnaseq_data)
corrected_sample_name <- substring(colnames(fibroblast_rnaseq_data),5,7)
colnames(fibroblast_rnaseq_data) <- corrected_sample_name

# Filter Lifespan Data to hFB12,13,14 
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_metadata <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_metadata <-  filtered_metadata[filtered_metadata$Treatment == "Galactose_21",]
filtered_metadata <-  filtered_metadata[filtered_metadata$Study_Part %in% c(4),]
filtered_metadata <-  filtered_metadata[!is.na(filtered_metadata$RNAseq_ID),]
filtered_metadata$RNAseq_ID
nrow(filtered_metadata)
# Retrieve TPM for control and Oligo samples
filtered_genes <- fibroblast_rnaseq_data[,colnames(fibroblast_rnaseq_data) %in% as.vector(filtered_metadata$RNAseq_ID)]
# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_metadata$RNAseq_ID, colnames(filtered_genes))]
#Timepoints
fibroblast_timepoints <- filtered_metadata$Days_Grown
# convert days grown to years
Timepoints_years <- fibroblast_timepoints /365
filtered_metadata <- data.frame(filtered_metadata, Timepoints_years)
## Check metadata matches data
all(colnames(filtered_genes)==filtered_metadata$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_metadata$RNAseq_ID))

galactose_fibroblasts_RNAseq_genes <- log2(filtered_genes+1)
galactose_fibroblasts_RNAseq_genes[1:5,1:5]
dim(galactose_fibroblasts_RNAseq_genes)
max(galactose_fibroblasts_RNAseq_genes)
min(galactose_fibroblasts_RNAseq_genes)
hist(as.matrix(galactose_fibroblasts_RNAseq_genes), main="Fibroblasts galactose RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)
galactose_fibroblasts_RNAseq_metadata <- filtered_metadata
setwd(baseDir)
save(galactose_fibroblasts_RNAseq_genes,file="Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_genes.RData")
write.csv(galactose_fibroblasts_RNAseq_metadata,"Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_metadata.csv")
#########

##### iPOP RNAseq ######
iPOP_visit_info <- read.csv("iPOP/iPOP_HMP&Exercise_visit_info_2.csv")
iPOP_study_info <- read.csv("iPOP/iPOP_HMP&Exercise_participant_info.csv")

iPOP_rnaseq_hmp <- read.csv("iPOP/iPOP_RNAseq_HMP_abundance_2.csv", row.names = 1)
iPOP_rnaseq_hmp[1:5,1:5]
dim(iPOP_rnaseq_hmp)
max(iPOP_rnaseq_hmp)
min(iPOP_rnaseq_hmp)
hist(as.matrix(iPOP_rnaseq_hmp), main="iPOP RNAseq Histogram", xlab = "log2 counts+1", breaks = 100)

# Get Age of subjects for each iPOP RNAseq sample
visitIDs <- gsub("-",".",iPOP_visit_info$VisitID)
rna_subjectIDs <- substring(colnames(iPOP_rnaseq_hmp), 1, str_locate(colnames(iPOP_rnaseq_hmp), "\\.")[,1]-1)
head(rna_subjectIDs)
length(rna_subjectIDs)
age_info <- data.frame(matrix(nrow = ncol(iPOP_rnaseq_hmp), ncol = 3))
colnames(age_info) <- c("visitID","subjectID","age")
for(i in 1:length(rna_subjectIDs)) {
  #print(i)
  subjectID <- rna_subjectIDs[i]
  #print(subjectID)
  age_info[i,2] <- subjectID

  visitID <- colnames(iPOP_rnaseq_hmp)[i]
  #print(visitID)
  age_info[i,1] <- visitID
  
  sampleAge <- iPOP_visit_info[visitIDs == visitID,]$combined_age
  #print(sampleAge)
  if(!visitID %in% visitIDs) {
    sampleAge <- iPOP_study_info[iPOP_study_info$SubjectID == subjectID,]$Age 
  }
  #print(sampleAge)
  age_info[i,3] <- sampleAge
  
}
head(age_info)
dim(age_info)

# Check that they match
all(colnames(iPOP_rnaseq_hmp)==age_info$visitID)
identical(colnames(iPOP_rnaseq_hmp),age_info$visitID)

setwd(baseDir)
iPOP_RNAseq_genes <- iPOP_rnaseq_hmp
iPOP_RNAseq_metadata <- age_info
save(iPOP_RNAseq_genes,file="iPOP/iPOP_RNAseq_genes.RData")
write.csv(iPOP_RNAseq_metadata,"iPOP/iPOP_RNAseq_metadata.csv")
#########

##### NIA RNAseq ######
AM_study_info <- read.csv("NIA/aging_muscle_metadata.csv")
AM_study_info$Sample
nSamples <- length(AM_study_info$Sample) # 53 samples
nSamples

AM_rnaseq_hmp <- read.csv("NIA/Aging_muscle_RNAseq_TPM.csv", row.names = 1)
AM_rnaseq_hmp <- log2(AM_rnaseq_hmp+1)
AM_rnaseq_hmp[1:5,1:5]
dim(AM_rnaseq_hmp)
max(AM_rnaseq_hmp)
min(AM_rnaseq_hmp)
hist(as.matrix(AM_rnaseq_hmp), main="NIA RNAseq Histogram", xlab = "log2 counts+1", breaks = 1000)

# Rename AM genes to gene symbols
fullgenenames <- rownames(AM_rnaseq_hmp)
genesymbols <- substring(fullgenenames,33,str_locate(fullgenenames, "protein_coding")[,1]-2)
#genesymbols <- substring(genesymbols, str_locate(genesymbols, "|")[,1])
head(genesymbols)
# Remove duplicate genes and genes with NA
AM_rnaseq_hmp <- AM_rnaseq_hmp[!duplicated(genesymbols, incomparables=""),]
dim(AM_rnaseq_hmp)
genesymbols <- genesymbols[!duplicated(genesymbols, incomparables="")]
genesymbols[is.na(genesymbols)] <- "klsdaj"
length(genesymbols)
rownames(AM_rnaseq_hmp) <- genesymbols

setwd(baseDir)
NIA_RNAseq_genes <- AM_rnaseq_hmp
NIA_RNAseq_metadata <- AM_study_info
save(NIA_RNAseq_genes,file="NIA/NIA_RNAseq_genes.RData")
write.csv(NIA_RNAseq_metadata,"NIA/NIA_RNAseq_metadata.csv")
#########

##### TwinsUK Skin RNAarray #######
# load in metadata
tUK_ages <- read.csv("TwinsUK/RNAarray/RNAarray_age_values.csv")
tUK_ages <- tUK_ages[order(tUK_ages$PUBLIC_ID, decreasing = F),]
dim(tUK_ages)

tUK_RNA_skin <- read.csv("TwinsUK/RNAarray/MuTHER_SKIN_main_normalized_data_newIds.csv", row.names = 1) 
tUK_RNA_skin[1:5,1:5]
max(tUK_RNA_skin)
min(tUK_RNA_skin)
dim(tUK_RNA_skin)

Gene_symbols <- data.frame(Gene=unlist(mget(x = rownames(tUK_RNA_skin),envir = illuminaHumanv3SYMBOL)))
head(Gene_symbols)
nrow(Gene_symbols)
length(unique(Gene_symbols[,1]))

# Filter for just unique genes
Gene <- Gene_symbols[,1]
tUK_RNA_skin_2 <- cbind(tUK_RNA_skin, Gene)
tUK_RNA_skin_2 <- tUK_RNA_skin_2[!is.na(tUK_RNA_skin_2$Gene),]
dim(tUK_RNA_skin_2)
genes <- tUK_RNA_skin_2$Gene
#rownames(tUK_RNA_skin_2) <- tUK_RNA_skin_2$Gene
tUK_RNA_skin_2 <- as.matrix(tUK_RNA_skin_2[,-706])

tUK_RNA_skin_3 <- matrix(ncol=ncol(tUK_RNA_skin_2),nrow=length(unique(genes)))
for(i in 1:length(unique(genes))) {
  gene <- unique(genes)[i]
  transcripts <- tUK_RNA_skin_2[genes == gene,]
  sum_gene <- transcripts
  if(is.matrix(transcripts)) {
    sum_gene <- colSums(transcripts)
  }
  tUK_RNA_skin_3[i,] <- sum_gene
  #rownames(tUK_RNA_skin_3)[i] <- ungene
  print(i)
}
colnames(tUK_RNA_skin_3) <- colnames(tUK_RNA_skin_2)
rownames(tUK_RNA_skin_3) <- unique(genes)
tUK_RNA_skin_3[1:5,1:5]
max(tUK_RNA_skin_3)
min(tUK_RNA_skin_3)
dim(tUK_RNA_skin_3)
hist(as.matrix(tUK_RNA_skin_3), main="TwinsUK SKin RNAarray Histogram", xlab = "expression", breaks = 100)

# order metadata and gene values
sample_donor <- substring(colnames(tUK_RNA_skin_3),1, str_locate(colnames(tUK_RNA_skin_3), "_")[,1]-1)
tUK_ages <- tUK_ages[tUK_ages$PUBLIC_ID %in% sample_donor,]
tUK_ages <- tUK_ages[match(sample_donor, tUK_ages$PUBLIC_ID),]
colnames(tUK_RNA_skin_3) <- sample_donor
identical(colnames(tUK_RNA_skin_3),tUK_ages$PUBLIC_ID)

setwd(baseDir)
twinsUK_skin_RNAarray_genes <- tUK_RNA_skin_3
twinsUK_skin_RNAarray_metadata <- tUK_ages
save(twinsUK_skin_RNAarray_genes,file="TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData")
write.csv(twinsUK_skin_RNAarray_metadata,"TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv")
#######


##### TwinsUK Fat RNAarray #######
# load in metadata
tUK_ages <- read.csv("TwinsUK/RNAarray/RNAarray_age_values.csv")
tUK_ages <- tUK_ages[order(tUK_ages$PUBLIC_ID, decreasing = F),]
dim(tUK_ages)

tUK_RNA_fat <- read.csv("TwinsUK/RNAarray/MuTHER_Fat_normalized_31032010_uncondensed_Ids.csv", row.names = 1) 
tUK_RNA_fat[1:5,1:5]
max(tUK_RNA_fat)
min(tUK_RNA_fat)
dim(tUK_RNA_fat)

Gene_symbols <- data.frame(Gene=unlist(mget(x = rownames(tUK_RNA_fat),envir = illuminaHumanv3SYMBOL)))
head(Gene_symbols)
nrow(Gene_symbols)
length(unique(Gene_symbols[,1]))

# Filter for just unique genes
Gene <- Gene_symbols[,1]
tUK_RNA_fat_2 <- cbind(tUK_RNA_fat, Gene)
tUK_RNA_fat_2 <- tUK_RNA_fat_2[!is.na(tUK_RNA_fat_2$Gene),]
dim(tUK_RNA_fat_2)
genes <- tUK_RNA_fat_2$Gene
#rownames(tUK_RNA_fat_2) <- tUK_RNA_fat_2$Gene
tUK_RNA_fat_2 <- as.matrix(tUK_RNA_fat_2[,-826])

tUK_RNA_fat_3 <- matrix(ncol=ncol(tUK_RNA_fat_2),nrow=length(unique(genes)))
for(i in 1:length(unique(genes))) {
  gene <- unique(genes)[i]
  transcripts <- tUK_RNA_fat_2[genes == gene,]
  sum_gene <- transcripts
  if(is.matrix(transcripts)) {
    sum_gene <- colSums(transcripts)
  }
  tUK_RNA_fat_3[i,]  <- sum_gene
  rownames(tUK_RNA_fat_3)[i] <- unique(genes)[i]
  #print(i)
}
colnames(tUK_RNA_fat_3) <- colnames(tUK_RNA_fat_2)
rownames(tUK_RNA_fat_3) <- unique(genes)
#tUK_RNA_fat_3 <- tUK_RNA_fat_2[!duplicated(tUK_RNA_fat_2$Gene),]
tUK_RNA_fat_3[1:5,1:5]
max(tUK_RNA_fat_3)
min(tUK_RNA_fat_3)
dim(tUK_RNA_fat_3)
hist(as.matrix(tUK_RNA_fat_3), main="TwinsUK Fat RNAarray Histogram", xlab = "expression", breaks = 100)

# order metadata and gene values
sample_donor <- substring(colnames(tUK_RNA_fat_3),1, str_locate(colnames(tUK_RNA_fat_3), "_")[,1]-1)
tUK_ages <- tUK_ages[tUK_ages$PUBLIC_ID %in% sample_donor,]
tUK_ages <- tUK_ages[match(sample_donor, tUK_ages$PUBLIC_ID),]
colnames(tUK_RNA_fat_3) <- sample_donor
identical(colnames(tUK_RNA_fat_3),tUK_ages$PUBLIC_ID)

setwd(baseDir)
twinsUK_fat_RNAarray_genes <- tUK_RNA_fat_3
twinsUK_fat_RNAarray_metadata <- tUK_ages
save(twinsUK_fat_RNAarray_genes,file="TwinsUK/RNAarray/twinsUK_fat_RNAarray_genes.RData")
write.csv(twinsUK_fat_RNAarray_metadata,"TwinsUK/RNAarray/twinsUK_fat_RNAarray_metadata.csv")
#######

##### Tabula Muris All Tissues RNAarray #######
# load in metadata
Tabula_Muris_metadata <- read.csv("Tabula_Muris/Tabula_Muris_metadata.csv",row.names=1)
head(Tabula_Muris_metadata)
dim(Tabula_Muris_metadata)
# load in data
Tabula_Muris_RNAseq_genes <- read.csv("Tabula_Muris/Tabula_Muris_bulkRNAseq_counttable_raw.csv",row.names=1)
Tabula_Muris_RNAseq_genes <- log2(Tabula_Muris_RNAseq_genes+1)
Tabula_Muris_RNAseq_genes[1:5,1:5]
max(Tabula_Muris_RNAseq_genes)
min(Tabula_Muris_RNAseq_genes)
dim(Tabula_Muris_RNAseq_genes)
hist(as.matrix(Tabula_Muris_RNAseq_genes), main="Tabula Muris RNAseq Histogram", xlab = "expression", breaks = 100)

# library(EnsDb.Mmusculus.v79)
# # 1. Convert from ensembl.gene to gene.symbol
# ensembl.genes <- rownames(Tabula_Muris_RNAseq_genes)
# geneIDs <- ensembldb::select(EnsDb.Mmusculus.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
# head(geneIDs)
# dim(geneIDs)

library(org.Mm.eg.db)
ensembl.genes <- rownames(Tabula_Muris_RNAseq_genes)
mapGenes <- mapIds(org.Mm.eg.db, ensembl.genes, keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
length(mapGenes)
identical(names(mapGenes), rownames(Tabula_Muris_RNAseq_genes))

# Filter for just unique genes
Tabula_Muris_RNAseq_genes <- cbind(Tabula_Muris_RNAseq_genes, mapGenes)
Tabula_Muris_RNAseq_genes <- Tabula_Muris_RNAseq_genes[!is.na(Tabula_Muris_RNAseq_genes$mapGenes),]
dim(Tabula_Muris_RNAseq_genes)
genes <- Tabula_Muris_RNAseq_genes$mapGenes
Tabula_Muris_RNAseq_genes <- as.matrix(Tabula_Muris_RNAseq_genes[,-893])

Tabula_Muris_RNAseq_genes_2 <- matrix(ncol=ncol(Tabula_Muris_RNAseq_genes),nrow=length(unique(genes)))
for(i in 1:length(unique(genes))) {
  gene <- unique(genes)[i]
  transcripts <- Tabula_Muris_RNAseq_genes[genes == gene,]
  sum_gene <- transcripts
  if(is.matrix(transcripts)) {
    sum_gene <- colSums(transcripts)
  }
  Tabula_Muris_RNAseq_genes_2[i,] <- sum_gene
  #rownames(Tabula_Muris_RNAseq_genes_2)[i] <- gene
  print(i)
}
colnames(Tabula_Muris_RNAseq_genes_2) <- colnames(Tabula_Muris_RNAseq_genes)
rownames(Tabula_Muris_RNAseq_genes_2) <- toupper(unique(genes))
Tabula_Muris_RNAseq_genes_2[1:5,1:5]
max(Tabula_Muris_RNAseq_genes_2)
min(Tabula_Muris_RNAseq_genes_2)
dim(Tabula_Muris_RNAseq_genes_2)
hist(as.matrix(Tabula_Muris_RNAseq_genes_2), main="Tabula Muris MultiTissue RNAseq Histogram", xlab = "expression", breaks = 100)


identical(gsub(".", "-", colnames(Tabula_Muris_RNAseq_genes), fixed=TRUE),as.character(Tabula_Muris_metadata$sample_id))


setwd(baseDir)
Tabula_Muris_RNAseq_genes <- Tabula_Muris_RNAseq_genes_2
save(Tabula_Muris_RNAseq_genes,file="Tabula_Muris/Tabula_Muris_RNAseq_genes.RData")
write.csv(Tabula_Muris_metadata,"Tabula_Muris/Tabula_Muris_RNAseq_metadata.csv")
#######

##### GTEx RNAseq #######
# load in metadata
GTEx_metadata <- read.csv("GTEx/GTEx_metadata.csv",row.names=1)
head(GTEx_metadata)
dim(GTEx_metadata)
# load in samples info
GTEx_sample_info <- read.csv("GTEx/GTEx_sample_info.csv",row.names=1)
head(GTEx_sample_info)
dim(GTEx_sample_info)
table(GTEx_sample_info$Note)
# load in data
#library(Rgraphviz) # BiocManager::install("Rgraphviz")
#library(CePa)
#GTEx_RNAseq_data <- read.gct("GTEx/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_rpkm.gct")
GTEx_RNAseq_data <- read.table("GTEx/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_rpkm.txt")
GTEx_RNAseq_data[1:5,1:5]
GTEx_RNAseq_transcript_names <- GTEx_RNAseq_data[,1]
GTEx_RNAseq_gene_names <- GTEx_RNAseq_data[,2]
length(unique(GTEx_RNAseq_gene_names))
rownames(GTEx_RNAseq_data) <-  GTEx_RNAseq_data[,1]
GTEx_RNAseq_data <- GTEx_RNAseq_data[,-1]
GTEx_RNAseq_data <- GTEx_RNAseq_data[,-1]
colnames(GTEx_RNAseq_data) <- GTEx_RNAseq_data[1,]
GTEx_RNAseq_data <- GTEx_RNAseq_data[-1,]
GTEx_RNAseq_data[1:5,1:5]

#write.csv(GTEx_RNAseq_data,"GTEx/GTEx_RNAseq_allTissues_gene_rpkm.csv")
GTEx_RNAseq_names <- data.frame(GTEx_RNAseq_transcript_names,GTEx_RNAseq_gene_names)
write.csv(GTEx_RNAseq_names,"GTEx/GTEx_RNAseq_transcript&gene_names.csv")

save(GTEx_RNAseq_matrix,file="GTEx/GTEx_RNAseq_allTissues_gene_rpkm.RData")

sum(as.numeric(GTEx_RNAseq_data[1,]))
GTEx_RNAseq_matrix <- as.matrix(sapply(GTEx_RNAseq_data, as.numeric))
typeof(GTEx_RNAseq_matrix[1,])
rownames(GTEx_RNAseq_matrix) <- GTEx_RNAseq_transcript_names[2:length(GTEx_RNAseq_transcript_names)]
GTEx_RNAseq_matrix[1:5,1:5]

# Sum transcripts that map to same gene
genes <- GTEx_RNAseq_gene_names[2:length(GTEx_RNAseq_gene_names)]
GTEx_RNAseq_matrix_2 <- matrix(ncol=ncol(GTEx_RNAseq_matrix),nrow=length(unique(genes)))
for(i in 1:length(unique(genes))) {
  gene <- unique(genes)[i]
  transcripts <- GTEx_RNAseq_matrix[genes == gene,]
  sum_gene <- transcripts
  if(is.matrix(transcripts)) {
    sum_gene <- colSums(transcripts)
  }
  GTEx_RNAseq_matrix_2[i,] <- sum_gene
  print(i)
}
colnames(GTEx_RNAseq_matrix_2) <- colnames(GTEx_RNAseq_matrix)
rownames(GTEx_RNAseq_matrix_2) <- unique(genes)
GTEx_RNAseq_matrix_2[1:5,1:5]
max(GTEx_RNAseq_matrix_2)
min(GTEx_RNAseq_matrix_2)
dim(GTEx_RNAseq_matrix_2)

GTEx_RNAseq_matrix_2 <- log2(GTEx_RNAseq_matrix_2+1)
GTEx_RNAseq_matrix_2[1:5,1:5]
max(GTEx_RNAseq_matrix_2)
min(GTEx_RNAseq_matrix_2)
dim(GTEx_RNAseq_matrix_2)
#hist(as.matrix(GTEx_RNAseq_matrix_2), main="GTEx RNAseq Histogram", xlab = "expression", breaks = 100)
#hist(as.matrix(GTEx_RNAseq_matrix_2), main="GTEx MultiTissue RNAseq Histogram", xlab = "expression", breaks = 100)


metadata <- as.data.frame(matrix(ncol=ncol(GTEx_metadata),nrow=nrow(GTEx_sample_info)))
leftS <-  substring(rownames(GTEx_sample_info),str_locate(rownames(GTEx_sample_info), "-")[,1]+1)
subjects <- paste0(substring(rownames(GTEx_sample_info),0,str_locate(rownames(GTEx_sample_info), "-")[,1]),substring(leftS,0,str_locate(leftS, "-")[,1]-1))
for(i in 1:length(subjects)) {
  subject <- subjects[i]
  metadata[i,] <- GTEx_metadata[GTEx_metadata$SUBJID ==subject,]
}
colnames(metadata) <- colnames(GTEx_metadata)
ageSubjects <- data.frame(subjects,metadata)
GTEx_sample_info <- data.frame(GTEx_sample_info, ageSubjects)
GTEx_sample_info <- GTEx_sample_info[match(colnames(GTEx_RNAseq_matrix_2),rownames(GTEx_sample_info)),]
identical(colnames(GTEx_RNAseq_matrix_2),rownames(GTEx_sample_info))

setwd(baseDir)
GTEx_RNAseq_data <- GTEx_RNAseq_matrix_2
GTEx_metadata <- GTEx_sample_info
save(GTEx_RNAseq_data,file="GTEx/GTEx_RNAseq_data.RData")
write.csv(GTEx_metadata,"GTEx/GTEx_RNAseq_metadata.csv")
#######

######## PART 1 - MODELING ########

# LMER Function
lmer_function <- function(data, timepoints, random_effect, filename) {

  numCores <- detectCores() - 2
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  inputs <- 1:nrow(data)
  LMER_modeling <- function(i) {
    
    molecule <- rownames(data)[i]
    rowData <- as.numeric(data[i,])
    
    #Test Gene
    #rowData <- data[59,]
    
    pValue <- 1
    intercept <- 0
    slope <- 0
    r2 <- 0
    if(sd(rowData) > 0) {
      model_data <- data.frame(rowData, timepoints, random_effect)
      
      MEmodel = lmer(rowData ~ timepoints + (1+timepoints|random_effect),
                     data = model_data,  REML = F)
      summary <- summary(MEmodel)
      intercept <- coef(summary)[,"Estimate"][1]
      slope <- coef(summary)[,"Estimate"][2]
      
      # model_data %>% 
      #   # save predicted values
      #   mutate(pred_dist = fitted(MEmodel)) %>% 
      #   # graph
      #   ggplot(aes(x=fibroblast_timepoints, y=pred_dist, group=random_effect, color=random_effect)) + theme_classic() +
      #   geom_line(size=1) 
      
      null_model <- lmer(rowData ~ (1|random_effect),
                         data = model_data,  REML = F)
      # model_data %>% 
      #   # save predicted values
      #   mutate(pred_dist = fitted(null_model)) %>% 
      #   # graph
      #   ggplot(aes(x=fibroblast_timepoints, y=pred_dist, group=random_effect, color=random_effect)) + theme_classic() +
      #   geom_line(size=1) 
      anova <- anova(null_model, MEmodel)
      pValue <- anova$`Pr(>Chisq)`[2]
      r2 <- r.squaredGLMM(MEmodel)[2]
    }
    message(i)
    
    #Organize matrix
    results <- data.frame(molecule, pValue, slope, intercept, r2)
    return(results)
    
  }
    
  
  lmer_results <- foreach(i=inputs, .packages=c('lme4','MuMIn')) %dopar% {
    LMER_modeling(i)
  }
  
  on.exit(stopCluster(cl))
  
  lmer_results <- do.call(rbind, lmer_results)
  head(lmer_results)

  # Correct for multipule Comparisons
  fdr <- as.matrix(p.adjust(lmer_results$pValue, method = "fdr", n = length(lmer_results$pValue)))
  lmer_results <- data.frame(lmer_results, fdr)

  write.csv(lmer_results, filename)
  return(lmer_results)
}

# GTEx LMER Function
# https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet
GTEx_lmer_function <- function(data, timepoints, random_effects, filename) {
  
  numCores <- detectCores() - 2
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  inputs <- 1:nrow(data)
  LMER_modeling <- function(i) {
    
    molecule <- rownames(data)[i]
    rowData <- as.numeric(data[i,])
    
    #Test Gene
    #tissue_data,timepoints,random_effects,filename
    #rowData <- tissue_data[59,]
    
    pValue <- 1
    intercept <- 0
    slope <- 0
    r2 <- 0
    if(sd(rowData) > 0) {
      model_data <- data.frame(rowData, timepoints, random_effects)
      
      lmer_code <- "rowData ~ timepoints "
      for(i in 1:ncol(random_effects)) {
        random_effect_name <- colnames(random_effects)[i]
        chr = "+ ("
        num = "1"
        if(length(unique(random_effects[,i])) > 1) {
          if(is.numeric(random_effects[,i])) { num="0"}
          lmer_code <- paste0(lmer_code, paste0(chr,num,"+timepoints|",random_effect_name,")"))
        }
      }
      lmer_code
      length(rowData)
      MEmodel = lmer(lmer_code,
                     data = model_data,  REML = F)
      summary <- summary(MEmodel)
      intercept <- coef(summary)[,"Estimate"][1]
      slope <- coef(summary)[,"Estimate"][2]
      
      # model_data %>%
      #   # save predicted values
      #   mutate(pred_dist = fitted(MEmodel)) %>%
      #   # graph
      #   ggplot(aes(x=timepoints, y=pred_dist, group=random_effects[,1], color=random_effects[,1])) + theme_classic() +
      #   geom_line(size=1)
      # 
      null_code <- "rowData ~ "
      for(i in 1:ncol(random_effects)) {
        random_effect_name <- colnames(random_effects)[i]
        chr = "+ ("
        if(i ==1) {chr="("}
        num = "1"
        if(length(unique(random_effects[,i])) > 1) {
          if(is.numeric(random_effects[,i])) { num="0"}
          null_code <- paste0(null_code, paste0(chr,num,"|",random_effect_name,")"))
        }
      }
      null_code  
      null_model <- lmer(null_code,
                         data = model_data,  REML = F)
      # model_data %>% 
      #   # save predicted values
      #   mutate(pred_dist = fitted(null_model)) %>% 
      #   # graph
      #   ggplot(aes(x=fibroblast_timepoints, y=pred_dist, group=random_effect, color=random_effect)) + theme_classic() +
      #   geom_line(size=1) 
      anova <- anova(null_model, MEmodel)
      pValue <- anova$`Pr(>Chisq)`[2]
      r2 <- r.squaredGLMM(MEmodel)[2]
    }
    message(i)
    
    #Organize matrix
    results <- data.frame(molecule, pValue, slope, intercept, r2)
    return(results)
    
  }
  
  
  lmer_results <- foreach(i=inputs, .packages=c('lme4','MuMIn')) %dopar% {
    LMER_modeling(i)
  }
  
  on.exit(stopCluster(cl))
  
  lmer_results <- do.call(rbind, lmer_results)
  head(lmer_results)
  
  # Correct for multipule Comparisons
  fdr <- as.matrix(p.adjust(lmer_results$pValue, method = "fdr", n = length(lmer_results$pValue)))
  lmer_results <- data.frame(lmer_results, fdr)
  
  write.csv(lmer_results, filename)
  return(lmer_results)
}


# Linear Function
linear_function <- function(data, timepoints, filename) {
  
  numCores <- detectCores() - 2
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  inputs <- 1:nrow(data)
  linear_modeling <- function(i) {
    
    molecule <- rownames(data)[i]
    rowData <- as.numeric(data[i,])
    #Test Gene
    #rowData <- data[59,]
    
    modelData <- data.frame(rowData,timepoints)
    pValue <- 1
    intercept <- 0
    slope <- 0
    r2 <- 0
    if(sd(rowData) > 0) {
      model <- lm(rowData ~ timepoints, data= modelData)
      intercept <- coef(model)[1]
      slope <- coef(model)[2]
      sum <- summary(model)
      r2 <- sum$r.squared
      pValue <- sum$coefficients[2,4]
    }
    message(i)
    
    #Organize matrix
    results <- data.frame(molecule, pValue, slope, intercept, r2)
    return(results)
    
  }
  
  
  linear_results <- foreach(i=inputs, .packages=c('lme4','MuMIn')) %dopar% {
    linear_modeling(i)
  }
  
  on.exit(stopCluster(cl))
  
  linear_results <- do.call(rbind, linear_results)
  head(linear_results)
  
  # Correct for multipule Comparisons
  fdr <- as.matrix(p.adjust(linear_results$pValue, method = "fdr", n = length(linear_results$pValue)))
  linear_results <- data.frame(linear_results, fdr)
  
  write.csv(linear_results, filename)
  return(linear_results)
}

# Permutation Function
set.seed(10)
nBatchSampling <- 50
nBootstrap <- 1000
percentile <- 0.01
Permutation_function <- function(data, timepoints, filename,nBoostrap,nBatchSampling,percentile) {
  # sampling_percentage = 0.6667
  #nSamples <- length(timepoints) * sampling_percentage
  numCores <- detectCores() - 2
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  inputs <- 1:nrow(data)
  Permuation_modeling <- function(i) {
    molecule <- rownames(data)[i]
    rowData <- as.numeric(data[i,])
    #Test Gene
    #rowData <- data[59,]
    # random sample subset of data
    modelData <- data.frame(rowData,timepoints)
    modelData <- modelData[sample(length(rowData), nBatchSampling),]
    boostrap_results <- matrix(ncol=4,nrow=nBoostrap)
    for(i in 1:nBoostrap) {
      pValue <- 1
      intercept <- 0
      slope <- 0
      r2 <- 0
      if(sd(modelData[,1]) > 0) {
        model <- lm(modelData[,1] ~ modelData[,2], data= modelData)
        intercept <- coef(model)[1]
        slope <- coef(model)[2]
        sum <- summary(model)
        r2 <- sum$r.squared
        pValue <- sum$coefficients[2,4]
      }
      #Organize matrix
      boostrap_results[i,] <- c(pValue, slope, intercept, r2)
    }
    # get 95% results
    #colnames(boostrap_results) <- c('pvalue','slope','intercept','r2')
    boostrap_results <- boostrap_results[order(abs(boostrap_results[,2]),decreasing=T),]
    results <- colMedians(boostrap_results[1:nBoostrap*percentile,])
    results <- c(molecule, results)
    return(results)
  }
  
  
  permutation_results <- foreach(i=inputs, .packages=c('lme4','MuMIn','robustbase')) %dopar% {
    Permuation_modeling(i)
  }
  
  on.exit(stopCluster(cl))
  
  permutation_results <- do.call(rbind, permutation_results)
  #permutation_results <- as.data.frame(permutation_results)
  colnames(permutation_results) <- c('molecule','pValue','slope','intercept','r2')
  #permutation_results$pValue <- as.numeric(permutation_results$pValue)
  #permutation_results$slope <- as.numeric(permutation_results$slope)
  #permutation_results$intercept <- as.numeric(permutation_results$intercept)
  #permutation_results$r2 <- as.numeric(permutation_results$r2)
  #head(permutation_results)
  
  # Correct for multipule Comparisons
  fdr <- p.adjust(as.numeric(permutation_results[,2]), method = "fdr", n = nrow(permutation_results)) # vector(mode='numeric',length=nrow(permutation_results))
  permutation_results <- data.frame(permutation_results, fdr)
  
  write.csv(permutation_results, filename)
  return(permutation_results)
}

# MinMax Function
MinMax_function <- function(data, timepoints, filename) {
  numCores <- detectCores() - 2
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  
  inputs <- 1:nrow(data)
  MinMax_modeling <- function(i) {
    
    molecule <- rownames(data)[i]
    rowData <- as.numeric(data[i,])
    #Test Gene
    #rowData <- data[59,]
    
    modelData <- data.frame(rowData,timepoints)
    pValue <- 1
    intercept <- 0
    slope <- 0
    r2 <- 0
    if(sd(rowData) > 0) {
      #min <- min(rowData)
      #minTimepoint <- timepoints[which.min(rowData)]
      #max <- max(rowData)
      #maxTimepoint <- timepoints[which.max(rowData)]
      #if(maxTimepoint < minTimepoint) {}
      #slope <- (max - min) / (maxTimepoint - minTimepoint)
      #intercept <- min -  slope*minTimepoint
      model <- lm(rowData ~ timepoints, data= modelData)
      intercept <- coef(model)[1]
      #slope <- coef(model)[2]
      slopes <- c(confint(model)[2,1],confint(model)[2,2])
      slope <- slopes[which(abs(slopes)==max(abs(slopes)))]
      sum <- summary(model)
      r2 <- sum$r.squared
      pValue <- sum$coefficients[2,4]
    }
    message(i)
    
    #Organize matrix
    results <- data.frame(molecule, pValue, slope, intercept, r2)
    return(results)
    
  }
  
  
  MinMax_results <- foreach(i=inputs, .packages=c('lme4','MuMIn')) %dopar% {
    MinMax_modeling(i)
  }
  
  on.exit(stopCluster(cl))
  
  MinMax_results <- do.call(rbind, MinMax_results)
  head(MinMax_results)
  
  # Correct for multipule Comparisons
  fdr <- as.matrix(p.adjust(MinMax_results$pValue, method = "fdr", n = length(MinMax_results$pValue)))
  MinMax_results <- data.frame(MinMax_results, fdr)
  write.csv(MinMax_results, filename)
  return(MinMax_results)
}


########### Linear Model for Control Fibroblasts r1 DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_betas.RData")
control_fibroblasts_r1_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_betas,timepoints,filename)
head(results)
dim(results)
###############

########### Linear Model for Control Fibroblasts r1-linear DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_linear_betas.RData")
control_fibroblasts_r1_linear_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_linear_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_linear_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_linear_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_linear_betas,timepoints,filename)
head(results)
dim(results)
###############

########### Linear Model for Control Fibroblasts r1-senescent DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescent_betas.RData")
control_fibroblasts_r1_senescent_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescent_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_senescent_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_senescent_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_senescent_betas,timepoints,filename)
head(results)
dim(results)
###############


########### Linear Model for Control Fibroblasts r1-earlyPhase DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_betas.RData")
control_fibroblasts_r1_earlyPhase_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_earlyPhase_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_earlyPhase_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_earlyPhase_betas,timepoints,filename)
head(results)
dim(results)
###############

########### Linear Model for Control Fibroblasts r1-midPhase DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_betas.RData")
control_fibroblasts_r1_midPhase_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_midPhase_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_midPhase_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_midPhase_betas,timepoints,filename)
head(results)
dim(results)
###############

########### Linear Model for Control Fibroblasts r1-senescentPhase DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_betas.RData")
control_fibroblasts_r1_senescentPhase_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_senescentPhase_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_senescentPhase_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_senescentPhase_betas,timepoints,filename)
head(results)
dim(results)
###############

########### Linear Model for Control Fibroblasts r1-treatments DNAm ######
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_betas.RData")
control_fibroblasts_r1_earlyPhase_betas[1:5,1:5]
control_fibroblasts_r1_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_metadata.csv")
control_fibroblasts_r1_metadata[1:5,1:5]
nrow(control_fibroblasts_r1_metadata)
identical(colnames(control_fibroblasts_r1_earlyPhase_betas),as.character(control_fibroblasts_r1_metadata$basename))
timepoints <- control_fibroblasts_r1_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r1_treatments_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_earlyPhase_betas,timepoints,filename)
head(results)
dim(results)
###############

# Linear Model for Control Fibroblasts r2 DNAm ###
##########
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_betas.RData")
control_fibroblasts_r2_betas[1:5,1:5]
control_fibroblasts_r2_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_metadata.csv")
nrow(control_fibroblasts_r2_metadata)
control_fibroblasts_r2_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r2_betas),as.character(control_fibroblasts_r2_metadata$basename))
timepoints <- control_fibroblasts_r2_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r2_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r2_linear_results.csv"
results <- linear_function(control_fibroblasts_r2_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Control Fibroblasts r3 DNAm ###
##########
load("Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_betas.RData")
control_fibroblasts_r3_betas[1:5,1:5]
control_fibroblasts_r3_metadata <- read.csv("Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_metadata.csv")
nrow(control_fibroblasts_r3_metadata)
control_fibroblasts_r3_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r3_betas),as.character(control_fibroblasts_r3_metadata$basename))
timepoints <- control_fibroblasts_r3_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r3_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/control_fibroblasts_r3_linear_results.csv"
results <- linear_function(control_fibroblasts_r3_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Contact Inhibited Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_betas.RData")
Contact_Inhibited_fibroblasts_betas <- filtered_fbr_betas
Contact_Inhibited_fibroblasts_betas[1:5,1:5]
Contact_Inhibited_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_metadata.csv")
Contact_Inhibited_fibroblasts_metadata[1:5,1:5]
identical(colnames(Contact_Inhibited_fibroblasts_betas),as.character(Contact_Inhibited_fibroblasts_metadata$basename))
timepoints <- Contact_Inhibited_fibroblasts_metadata$Timepoints_years
cell_lines <- Contact_Inhibited_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/Contact_Inhibition_fibroblasts_linear_results.csv"
Contact_Inhibited_fibroblasts_betas_lmer_results <- linear_function(Contact_Inhibited_fibroblasts_betas,timepoints,filename)
head(Contact_Inhibited_fibroblasts_betas_lmer_results)
dim(Contact_Inhibited_fibroblasts_betas_lmer_results)
##########

# Linear Model for Contact Inhibited earlyPhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_betas.RData")
Contact_Inhibited_fibroblasts_betas <- filtered_fbr_betas
Contact_Inhibited_fibroblasts_betas[1:5,1:5]
metadata <- read.csv("Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_metadata.csv")
metadata[1:5,1:5]
nrow(metadata)
identical(colnames(Contact_Inhibited_fibroblasts_betas),as.character(metadata$basename))
timepoints <- metadata$Timepoints_years
cell_lines <- metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_fibroblasts_linear_results.csv"
results <- linear_function(Contact_Inhibited_fibroblasts_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Contact Inhibited latePhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_betas.RData")
Contact_Inhibited_fibroblasts_betas <- filtered_fbr_betas
Contact_Inhibited_fibroblasts_betas[1:5,1:5]
metadata <- read.csv("Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_metadata.csv")
metadata[1:5,1:5]
nrow(metadata)
identical(colnames(Contact_Inhibited_fibroblasts_betas),as.character(metadata$basename))
timepoints <- metadata$Timepoints_years
cell_lines <- metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_fibroblasts_linear_results.csv"
results <- linear_function(Contact_Inhibited_fibroblasts_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Hypoxia Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_betas.RData")
hypoxia_betas[1:5,1:5]
Hypoxia_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_metadata.csv")
Hypoxia_fibroblasts_metadata[1:5,1:5]
identical(colnames(hypoxia_betas),as.character(Hypoxia_fibroblasts_metadata$basename))
timepoints <- Hypoxia_fibroblasts_metadata$Timepoints_years
cell_lines <- Hypoxia_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/Hypoxia_fibroblasts_linear_results.csv"
results <- linear_function(hypoxia_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for DEX earlyPhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_betas.RData")
DEX_betas[1:5,1:5]
DEX_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_metadata.csv")
DEX_fibroblasts_metadata[1:5,1:5]
nrow(DEX_fibroblasts_metadata)
identical(colnames(DEX_betas),as.character(DEX_fibroblasts_metadata$basename))
timepoints <- DEX_fibroblasts_metadata$Timepoints_years
cell_lines <- DEX_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/DEX_earlyPhase_fibroblasts_linear_results.csv"
results <- linear_function(DEX_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Modulators earlyPhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_betas.RData")
modulators_betas[1:5,1:5]
modulators_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_metadata.csv")
modulators_fibroblasts_metadata[1:5,1:5]
nrow(modulators_fibroblasts_metadata)
identical(colnames(modulators_betas),as.character(modulators_fibroblasts_metadata$basename))
timepoints <- modulators_fibroblasts_metadata$Timepoints_years
cell_lines <- modulators_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/modulators_earlyPhase_fibroblasts_linear_results.csv"
results <- linear_function(modulators_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for Oligomycin earlyPhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_betas.RData")
oligo_betas[1:5,1:5]
oligo_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_metadata.csv")
oligo_fibroblasts_metadata[1:5,1:5]
nrow(oligo_fibroblasts_metadata)
identical(colnames(oligo_betas),as.character(oligo_fibroblasts_metadata$basename))
timepoints <- oligo_fibroblasts_metadata$Timepoints_years
cell_lines <- oligo_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/oligo_earlyPhase_fibroblasts_linear_results.csv"
results <- linear_function(oligo_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for SURF1 earlyPhase Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_betas.RData")
SURF1_betas[1:5,1:5]
SURF1_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_metadata.csv")
SURF1_fibroblasts_metadata[1:5,1:5]
nrow(SURF1_fibroblasts_metadata)
identical(colnames(SURF1_betas),as.character(SURF1_fibroblasts_metadata$basename))
timepoints <- SURF1_fibroblasts_metadata$Timepoints_years
cell_lines <- SURF1_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/SURF1_earlyPhase_fibroblasts_linear_results.csv"
results <- linear_function(SURF1_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for 2-Deoxyglucose Fibroblasts DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_betas.RData")
deoxyglucose_betas[1:5,1:5]
deoxyglucose_fibroblasts_metadata <- read.csv("Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_metadata.csv")
deoxyglucose_fibroblasts_metadata[1:5,1:5]
identical(colnames(deoxyglucose_betas),as.character(deoxyglucose_fibroblasts_metadata$basename))
timepoints <- deoxyglucose_fibroblasts_metadata$Timepoints_years
cell_lines <- deoxyglucose_fibroblasts_metadata$Cell_Line
filename <- "Cellular_Lifespan/DNAm/deoxyglucose_fibroblasts_linear_results.csv"
results <- linear_function(deoxyglucose_betas,timepoints,filename)
head(results)
dim(results)
##########

# Linear Model for HEK293 DNAm
##########
setwd(baseDir)
load("Cellular_Lifespan/DNAm/HEK293_betas.RData")
HEK293_betas[1:5,1:5]
HEK293_metadata <- read.csv("Cellular_Lifespan/DNAm/HEK293_metadata.csv")
HEK293_metadata[1:5,1:5]
identical(colnames(HEK293_betas),as.character(HEK293_metadata$basename))
timepoints <- HEK293_metadata$Timepoints_years
filename <- "Cellular_Lifespan/DNAm/HEK293_DNAm_linear_results.csv"
HEK293_betas_linear_results <- linear_function(HEK293_betas,timepoints,filename)
head(HEK293_betas_linear_results)
dim(HEK293_betas_linear_results)
##########

# Linear  Model for TwinsUK Skin DNAm
##########
load("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData")
twinsUK_skin_DNAm_betas[1:5,1:5]
dim(twinsUK_skin_DNAm_betas)
twinsUK_skin_DNAm_metadata <- read.csv("TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv")
twinsUK_skin_DNAm_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_DNAm_betas),as.character(twinsUK_skin_DNAm_metadata$X))
timepoints <- twinsUK_skin_DNAm_metadata$age.at.biopsy.ch1
length(timepoints) == ncol(twinsUK_skin_DNAm_betas)
twins <- twinsUK_skin_DNAm_metadata$relatedness.identification.id.ch1
filename <- "TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_linear_results.csv"
twinsUK_skin_DNAm_linear_results <- linear_function(twinsUK_skin_DNAm_betas,timepoints,filename)
head(twinsUK_skin_DNAm_linear_results)
dim(twinsUK_skin_DNAm_linear_results)
##########

# Permutation  Model for TwinsUK Skin DNAm
##########
load("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData")
twinsUK_skin_DNAm_betas[1:5,1:5]
dim(twinsUK_skin_DNAm_betas)
twinsUK_skin_DNAm_metadata <- read.csv("TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv")
twinsUK_skin_DNAm_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_DNAm_betas),as.character(twinsUK_skin_DNAm_metadata$X))
timepoints <- twinsUK_skin_DNAm_metadata$age.at.biopsy.ch1
length(timepoints) == ncol(twinsUK_skin_DNAm_betas)
twins <- twinsUK_skin_DNAm_metadata$relatedness.identification.id.ch1
# limit to significant sites from linear model
linear_results <- read.csv('TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_linear_results.csv')
sig_genes <- linear_results[linear_results$fdr < 0.05,]$molecule
length(sig_genes)
twinsUK_skin_DNAm_betas <- twinsUK_skin_DNAm_betas[rownames(twinsUK_skin_DNAm_betas) %in% sig_genes,]
dim(twinsUK_skin_DNAm_betas)
print(identical(colnames(twinsUK_skin_DNAm_betas), as.character(twinsUK_skin_DNAm_metadata$X)))
filename <- "TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_permutation_results.csv"
results <- Permutation_function(twinsUK_skin_DNAm_betas,timepoints,filename,nBootstrap,nBatchSampling,percentile)
head(results)
dim(results)
##########

# MinMax  Model for TwinsUK Skin DNAm
##########
load("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData")
twinsUK_skin_DNAm_betas[1:5,1:5]
dim(twinsUK_skin_DNAm_betas)
twinsUK_skin_DNAm_metadata <- read.csv("TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv")
twinsUK_skin_DNAm_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_DNAm_betas),as.character(twinsUK_skin_DNAm_metadata$X))
timepoints <- twinsUK_skin_DNAm_metadata$age.at.biopsy.ch1
length(timepoints) == ncol(twinsUK_skin_DNAm_betas)
twins <- twinsUK_skin_DNAm_metadata$relatedness.identification.id.ch1
# limit to significant sites from linear model
linear_results <- read.csv('TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_linear_results.csv')
sig_genes <- linear_results[linear_results$fdr < 0.05,]$molecule
length(sig_genes)
twinsUK_skin_DNAm_betas <- twinsUK_skin_DNAm_betas[rownames(twinsUK_skin_DNAm_betas) %in% sig_genes,]
dim(twinsUK_skin_DNAm_betas)
print(identical(colnames(twinsUK_skin_DNAm_betas), as.character(twinsUK_skin_DNAm_metadata$X)))
filename <- "TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_minmax_results.csv"
results <- MinMax_function(twinsUK_skin_DNAm_betas,timepoints,filename)
head(results)
dim(results)
##########


# Linear  Model for TwinsUK Fat DNAm
##########
load("TwinsUK/DNAm_Fat/twinsUK_FAT_DNAm_betas.RData")
twinsUK_fat_DNAm_betas[1:5,1:5]
dim(twinsUK_fat_DNAm_betas)
twinsUK_Fat_DNAm_metadata <- read.csv("TwinsUK/DNAm_Fat/RNAarray_age_values.csv")
twinsUK_Fat_DNAm_metadata[1:4,1:4]
twinsUK_Fat_DNAm_metadata <- twinsUK_Fat_DNAm_metadata[twinsUK_Fat_DNAm_metadata$PUBLIC_ID %in% colnames(twinsUK_fat_DNAm_betas),]
twinsUK_Fat_DNAm_metadata <- twinsUK_Fat_DNAm_metadata[match(colnames(twinsUK_fat_DNAm_betas), twinsUK_Fat_DNAm_metadata$PUBLIC_ID),]
identical(colnames(twinsUK_fat_DNAm_betas),as.character(twinsUK_Fat_DNAm_metadata$PUBLIC_ID))
timepoints <- twinsUK_Fat_DNAm_metadata$AGE
length(timepoints) == ncol(twinsUK_fat_DNAm_betas)
twins <- twinsUK_Fat_DNAm_metadata$CO.TWIN
filename <- "TwinsUK/DNAm_Fat/twinsUK_Fat_DNAm_linear_results.csv"
twinsUK_Fat_DNAm_linear_results <- linear_function(twinsUK_fat_DNAm_betas,timepoints, filename)
head(twinsUK_Fat_DNAm_linear_results)
dim(twinsUK_Fat_DNAm_linear_results)
##########

############

## Linear Model for Control Fibroblasts r1 RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1-linear RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1-senescent RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1 earlyPhase RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1 midPhase RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midphase_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1 senescentPhase RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r1 treatments RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_genes.RData")
control_fibroblasts_r1_RNAseq_genes[1:5,1:5]
control_fibroblasts_r1_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_metadata.csv")
control_fibroblasts_r1_RNAseq_metadata[1:5,1:5]
nrow(control_fibroblasts_r1_RNAseq_metadata)
identical(colnames(control_fibroblasts_r1_RNAseq_genes),as.character(control_fibroblasts_r1_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r1_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r1_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r1_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r2 RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_genes.RData")
control_fibroblasts_r2_RNAseq_genes[1:5,1:5]
control_fibroblasts_r2_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_metadata.csv")
control_fibroblasts_r2_RNAseq_metadata[1:5,1:5]
identical(colnames(control_fibroblasts_r2_RNAseq_genes),as.character(control_fibroblasts_r2_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r2_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r2_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r2_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r3 RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_genes.RData")
control_fibroblasts_r3_RNAseq_genes[1:5,1:5]
control_fibroblasts_r3_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_metadata.csv")
control_fibroblasts_r3_RNAseq_metadata[1:5,1:5]
nrow(control_fibroblasts_r3_RNAseq_metadata)
identical(colnames(control_fibroblasts_r3_RNAseq_genes),as.character(control_fibroblasts_r3_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r3_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r3_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r3_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r3 earlyPhase RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase_RNAseq_genes.RData")
control_fibroblasts_r3_RNAseq_genes[1:5,1:5]
control_fibroblasts_r3_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase_RNAseq_metadata.csv")
control_fibroblasts_r3_RNAseq_metadata[1:5,1:5]
nrow(control_fibroblasts_r3_RNAseq_metadata)
identical(colnames(control_fibroblasts_r3_RNAseq_genes),as.character(control_fibroblasts_r3_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r3_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r3_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r3_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Control Fibroblasts r3 earlyPhase13 RNAseq
##########
load("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase13_RNAseq_genes.RData")
control_fibroblasts_r3_RNAseq_genes[1:5,1:5]
control_fibroblasts_r3_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase13_RNAseq_metadata.csv")
control_fibroblasts_r3_RNAseq_metadata[1:5,1:5]
nrow(control_fibroblasts_r3_RNAseq_metadata)
identical(colnames(control_fibroblasts_r3_RNAseq_genes),as.character(control_fibroblasts_r3_RNAseq_metadata$RNAseq_ID))
timepoints <- control_fibroblasts_r3_RNAseq_metadata$Timepoints_years
cell_lines <- control_fibroblasts_r3_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/control_fibroblasts_r3_earlyPhase13_RNAseq_linear_results.csv"
results <- linear_function(control_fibroblasts_r3_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Contact_Inhibition Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/contact_inhibition_fibroblasts_RNAseq_genes.RData")
contact_inhibition_fibroblasts_RNAseq_genes[1:5,1:5]
contact_inhibition_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_RNAseq_metadata.csv")
contact_inhibition_fibroblasts_RNAseq_metadata[1:5,1:5]
identical(colnames(contact_inhibition_fibroblasts_RNAseq_genes),as.character(contact_inhibition_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- contact_inhibition_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- contact_inhibition_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/contact_inhibition_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(contact_inhibition_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Contact_Inhibition earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes.RData")
contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes[1:5,1:5]
contact_inhibition_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_RNAseq_metadata.csv")
contact_inhibition_fibroblasts_RNAseq_metadata[1:5,1:5]
identical(colnames(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes),as.character(contact_inhibition_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- contact_inhibition_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- contact_inhibition_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Contact_Inhibition latePhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_fibroblasts_RNAseq_genes.RData")
contact_inhibition_latePhase_fibroblasts_RNAseq_genes[1:5,1:5]
contact_inhibition_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_RNAseq_metadata.csv")
contact_inhibition_fibroblasts_RNAseq_metadata[1:5,1:5]
identical(colnames(contact_inhibition_latePhase_fibroblasts_RNAseq_genes),as.character(contact_inhibition_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- contact_inhibition_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- contact_inhibition_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(contact_inhibition_latePhase_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Hypoxia Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_genes.RData")
hypoxia_fibroblasts_RNAseq_genes[1:5,1:5]
hypoxia_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_metadata.csv")
hypoxia_fibroblasts_RNAseq_metadata[1:5,1:5]
identical(colnames(hypoxia_fibroblasts_RNAseq_genes),as.character(hypoxia_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- hypoxia_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- hypoxia_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(hypoxia_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Hypoxia Contact Inhibition Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_genes.RData")
hypoxia_contact_inhibition_fibroblasts_RNAseq_genes[1:5,1:5]
hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata.csv")
hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata)
identical(colnames(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes),as.character(hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(hypoxia_contact_inhibition_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for DEX Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_genes.RData")
dex_fibroblasts_RNAseq_genes[1:5,1:5]
dex_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_metadata.csv")
dex_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(dex_fibroblasts_RNAseq_metadata)
identical(colnames(dex_fibroblasts_RNAseq_genes),as.character(dex_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- dex_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- dex_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(dex_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for DEX earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_genes.RData")
dex_fibroblasts_RNAseq_genes[1:5,1:5]
dex_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_metadata.csv")
dex_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(dex_fibroblasts_RNAseq_metadata)
identical(colnames(dex_fibroblasts_RNAseq_genes),as.character(dex_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- dex_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- dex_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(dex_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for SURF1 Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_genes.RData")
surf1_fibroblasts_RNAseq_genes[1:5,1:5]
surf1_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_metadata.csv")
surf1_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(surf1_fibroblasts_RNAseq_metadata)
identical(colnames(surf1_fibroblasts_RNAseq_genes),as.character(surf1_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- surf1_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- surf1_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(surf1_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for SURF1 earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_genes.RData")
surf1_fibroblasts_RNAseq_genes[1:5,1:5]
surf1_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_metadata.csv")
surf1_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(surf1_fibroblasts_RNAseq_metadata)
identical(colnames(surf1_fibroblasts_RNAseq_genes),as.character(surf1_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- surf1_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- surf1_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(surf1_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Oligomycin Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_genes.RData")
oligo_fibroblasts_RNAseq_genes[1:5,1:5]
oligo_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_metadata.csv")
oligo_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(oligo_fibroblasts_RNAseq_metadata)
identical(colnames(oligo_fibroblasts_RNAseq_genes),as.character(oligo_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- oligo_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- oligo_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(oligo_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Oligomycin earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_genes.RData")
oligomycin_fibroblasts_RNAseq_genes[1:5,1:5]
oligo_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_metadata.csv")
oligo_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(oligo_fibroblasts_RNAseq_metadata)
identical(colnames(oligomycin_fibroblasts_RNAseq_genes),as.character(oligo_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- oligo_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- oligo_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/oligo_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(oligomycin_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Mito-Modulators Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_genes.RData")
modulators_fibroblasts_RNAseq_genes[1:5,1:5]
modulators_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_metadata.csv")
modulators_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(modulators_fibroblasts_RNAseq_metadata)
identical(colnames(modulators_fibroblasts_RNAseq_genes),as.character(modulators_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- modulators_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- modulators_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(modulators_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for Mito-Modulators earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_genes.RData")
modulators_fibroblasts_RNAseq_genes[1:5,1:5]
modulators_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_metadata.csv")
modulators_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(modulators_fibroblasts_RNAseq_metadata)
identical(colnames(modulators_fibroblasts_RNAseq_genes),as.character(modulators_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- modulators_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- modulators_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(modulators_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for 2-deoxyglucose Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_genes.RData")
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
deoxyglucose_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_metadata.csv")
deoxyglucose_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(deoxyglucose_fibroblasts_RNAseq_metadata)
identical(colnames(deoxyglucose_fibroblasts_RNAseq_genes),as.character(deoxyglucose_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- deoxyglucose_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- deoxyglucose_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(deoxyglucose_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for 2-deoxyglucose earlyPhase Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_genes.RData")
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
deoxyglucose_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_metadata.csv")
deoxyglucose_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(deoxyglucose_fibroblasts_RNAseq_metadata)
identical(colnames(deoxyglucose_fibroblasts_RNAseq_genes),as.character(deoxyglucose_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- deoxyglucose_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- deoxyglucose_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(deoxyglucose_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for 2-deoxyglucose earlyPhase13 Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase13_fibroblasts_RNAseq_genes.RData")
deoxyglucose_fibroblasts_RNAseq_genes[1:5,1:5]
deoxyglucose_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase13_fibroblasts_RNAseq_metadata.csv")
deoxyglucose_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(deoxyglucose_fibroblasts_RNAseq_metadata)
identical(colnames(deoxyglucose_fibroblasts_RNAseq_genes),as.character(deoxyglucose_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- deoxyglucose_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- deoxyglucose_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase13_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(deoxyglucose_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########


## Linear Model for betahydroxybutyrate Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_genes.RData")
betahydroxybutyrate_fibroblasts_RNAseq_genes[1:5,1:5]
betahydroxybutyrate_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_metadata.csv")
betahydroxybutyrate_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(betahydroxybutyrate_fibroblasts_RNAseq_metadata)
identical(colnames(betahydroxybutyrate_fibroblasts_RNAseq_genes),as.character(betahydroxybutyrate_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- betahydroxybutyrate_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- betahydroxybutyrate_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(betahydroxybutyrate_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for galactose Fibroblasts RNAseq
##########
load("Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_genes.RData")
galactose_fibroblasts_RNAseq_genes[1:5,1:5]
galactose_fibroblasts_RNAseq_metadata <- read.csv("Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_metadata.csv")
galactose_fibroblasts_RNAseq_metadata[1:5,1:5]
nrow(galactose_fibroblasts_RNAseq_metadata)
identical(colnames(galactose_fibroblasts_RNAseq_genes),as.character(galactose_fibroblasts_RNAseq_metadata$RNAseq_ID))
timepoints <- galactose_fibroblasts_RNAseq_metadata$Timepoints_years
cell_lines <- galactose_fibroblasts_RNAseq_metadata$Cell_Line
filename <- "Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_linear_results.csv"
results <- linear_function(galactose_fibroblasts_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## LMER Model for iPOP RNAseq
##########
load("iPOP/iPOP_RNAseq_genes.RData")
iPOP_RNAseq_genes[1:5,1:5]
dim(iPOP_RNAseq_genes)
iPOP_RNAseq_metadata <- read.csv("iPOP/iPOP_RNAseq_metadata.csv")
iPOP_RNAseq_metadata[1:5,1:3]
identical(colnames(iPOP_RNAseq_genes),as.character(iPOP_RNAseq_metadata$visitID))
timepoints <- iPOP_RNAseq_metadata$age
subjects <- iPOP_RNAseq_metadata$subjectID
filename <- "iPOP/iPOP_RNAseq_lmer_results.csv"
iPOP_RNAseq_lmer_results <- lmer_function(iPOP_RNAseq_genes,timepoints,subjects,filename)
head(iPOP_RNAseq_lmer_results)
dim(iPOP_RNAseq_lmer_results)
##########

## Linear Model for NIA RNAseq
##########
load("NIA/NIA_RNAseq_genes.RData")
NIA_RNAseq_genes[1:5,1:5]
dim(NIA_RNAseq_genes)
NIA_RNAseq_metadata <- read.csv("NIA/NIA_RNAseq_metadata.csv")
NIA_RNAseq_metadata[1:5,1:5]
identical(colnames(NIA_RNAseq_genes),as.character(NIA_RNAseq_metadata$Sample))
timepoints <- NIA_RNAseq_metadata$Age
filename <- "NIA/NIA_RNAseq_linear_results.csv"
results <- linear_function(NIA_RNAseq_genes,timepoints,filename)
head(results)
dim(results)
##########

## Linear Model for TwinsUK Skin RNAarray
##########
load("TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData")
twinsUK_skin_RNAarray_genes[1:5,1:5]
dim(twinsUK_skin_RNAarray_genes)
twinsUK_skin_RNAarray_metadata <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv")
twinsUK_skin_RNAarray_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_RNAarray_genes),as.character(twinsUK_skin_RNAarray_metadata$PUBLIC_ID))
timepoints <- twinsUK_skin_RNAarray_metadata$AGE
twins <- twinsUK_skin_RNAarray_metadata$CO.TWIN
filename <- "TwinsUK/RNAarray/twinsUK_skin_RNAarray_linear_results.csv"
results <- linear_function(twinsUK_skin_RNAarray_genes,timepoints,filename)
head(results)
dim(results)
##########

## ConfidenceInterval Model for TwinsUK Skin RNAarray
##########
load("TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData")
twinsUK_skin_RNAarray_genes[1:5,1:5]
dim(twinsUK_skin_RNAarray_genes)
twinsUK_skin_RNAarray_metadata <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv")
twinsUK_skin_RNAarray_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_RNAarray_genes),as.character(twinsUK_skin_RNAarray_metadata$PUBLIC_ID))
timepoints <- twinsUK_skin_RNAarray_metadata$AGE
twins <- twinsUK_skin_RNAarray_metadata$CO.TWIN
filename <- "TwinsUK/RNAarray/twinsUK_skin_RNAarray_minmax_results.csv"
results <- MinMax_function(twinsUK_skin_RNAarray_genes,timepoints,filename)
results <- results[order(abs(results$slope),decreasing=T),]
head(results)
dim(results)
##########


## Permutation Model for TwinsUK Skin RNAarray
##########
load("TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData")
twinsUK_skin_RNAarray_genes[1:5,1:5]
dim(twinsUK_skin_RNAarray_genes)
twinsUK_skin_RNAarray_metadata <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv")
twinsUK_skin_RNAarray_metadata[1:5,1:5]
identical(colnames(twinsUK_skin_RNAarray_genes),as.character(twinsUK_skin_RNAarray_metadata$PUBLIC_ID))
timepoints <- twinsUK_skin_RNAarray_metadata$AGE
twins <- twinsUK_skin_RNAarray_metadata$CO.TWIN
filename <- "TwinsUK/RNAarray/twinsUK_skin_RNAarray_permutation_results.csv"
nBatchSampling <- round(nrow(twinsUK_skin_RNAarray_metadata) * 0.1,0)
results <- Permutation_function(twinsUK_skin_RNAarray_genes,timepoints,filename,nBootstrap,nBatchSampling,percentile)
results <- results[order(abs(as.numeric(results$slope)),decreasing=T),]
head(results)
dim(results)
nrow(results[results$fdr <0.05,])
##########

## Linear Model for TwinsUK Fat RNAarray
##########
load("TwinsUK/RNAarray/twinsUK_fat_RNAarray_genes.RData")
twinsUK_fat_RNAarray_genes[1:5,1:5]
dim(twinsUK_fat_RNAarray_genes)
twinsUK_fat_RNAarray_metadata <- read.csv("TwinsUK/RNAarray/twinsUK_fat_RNAarray_metadata.csv")
twinsUK_fat_RNAarray_metadata[1:5,1:5]
identical(colnames(twinsUK_fat_RNAarray_genes),as.character(twinsUK_fat_RNAarray_metadata$PUBLIC_ID))
timepoints <- twinsUK_fat_RNAarray_metadata$AGE
twins <- twinsUK_fat_RNAarray_metadata$CO.TWIN
filename <- "TwinsUK/RNAarray/twinsUK_fat_RNAarray_linear_results.csv"
results <- linear_function(twinsUK_fat_RNAarray_genes,timepoints,filename)
head(results)
dim(results)
##########

## LMER Model for Tabula Muris MultiTissue RNAseq
##########
load("Tabula_Muris/Tabula_Muris_RNAseq_genes.RData")
Tabula_Muris_RNAseq_genes[1:5,1:5]
dim(Tabula_Muris_RNAseq_genes)
metadata <- read.csv("Tabula_Muris/Tabula_Muris_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(gsub(".", "-", colnames(Tabula_Muris_RNAseq_genes), fixed=TRUE),as.character(metadata$sample_id))
timepoints <- metadata$age
tissue <- metadata$tissue
filename <- "Tabula_Muris/Tabula_Muris_RNAseq_lmer_results.csv"
results <- lmer_function(Tabula_Muris_RNAseq_genes,timepoints,tissue,filename)
head(results)
dim(results)
##########

# grouping skin GTEx data
skin_exp_meta <- metadata[metadata$Note == 'Skin - Sun Exposed (Lower leg)',]
nrow(skin_exp_meta)
skin_cov_meta <- metadata[metadata$Note == 'Skin - Not Sun Exposed (Suprapubic)',]
nrow(skin_cov_meta)
length(intersect(skin_exp_meta$SUBJID,skin_cov_meta$SUBJID))

## LMER Model for GTEx MultiTissue RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
subTissues <- metadata$Note
table(subTissues)

# Break up by tissue
tissues <- c("Adipose","Adrenal_Gland","Artery","Brain","Breast",
             "Colon","Esophagus","Fibroblasts","Liver",
             "Heart","Lung","Muscle","Nerve","Ovary","Pancreas","Pituitary","Prostate",
             "Skin","Small_Intestine","Spleen","Stomach","Testis","Thyroid","Uterus","Vagina","Whole_Blood")
length(tissues[tissues %in% unique(metadata$Tissue)])
nTissues <- length(tissues)
nTissues
fileData <- matrix(ncol=5,nrow=nTissues)
colnames(fileData) <- c("Tissue","nSamples","filename","nSubTissues","Random_Effect_Category")
for(i in 1:nTissues) {
  tissue <- tissues[i]
  print(tissue)
  tissue_metadata <- metadata[metadata$Tissue == tissue,]
  tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
  print(identical(colnames(tissue_data), rownames(tissue_metadata)))
  timepoints <- tissue_metadata$AGE
  random_effect <- tissue_metadata$COHORT
  random_effect_category <- "DonationType"
  print(unique(tissue_metadata$Note))
  if(length(unique(tissue_metadata$Note)) >1) {
    random_effect <- tissue_metadata$Note
    random_effect_category <- "SubTissue"
  }
  else if(length(unique(random_effect)) == 1) {
    random_effect <- tissue_metadata$SEX
    random_effect_category <- "Sex"
  }
  else if(length(unique(random_effect)) == 1) {
    random_effect <- tissue_metadata$RACE
    random_effect_category <- "Race"
  }
  print(random_effect_category)
  print("modeling..")
  filename <- paste0("GTEx/GTEx",tissue,"_RNAseq_lmer_results.csv")
  #results <- lmer_function(tissue_data,timepoints,random_effect,filename)
  #head(results)
  #dim(results)
  nSubTissues <- length(unique(tissue_metadata$Note))
  fileData[i,] <- c(tissue, nrow(tissue_metadata),filename,nSubTissues,random_effect_category)
}
write.csv(fileData,"GTEx/tissueFileMetadata.csv")

##########

## GTEx LMER Model for GTEx MultiTissue RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
subTissues <- metadata$Note
table(subTissues)

# Break up by tissue
tissues <- unique(subTissues)
nTissues <- length(tissues)
nTissues
fileData <- matrix(ncol=7,nrow=nTissues)
colnames(fileData) <- c("Tissue","nSamples","filename","nSubTissues","datafilename","metadatafilename","model")
for(i in 1:nTissues) {
  tissue <- tissues[i]
  
  print(tissue)
  tissue_metadata <- metadata[metadata$Note == tissue,]
  if(nrow(tissue_metadata) < 50 | tissue == "Cells - Leukemia cell line (CML)") {
    next
  }
  tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
  # limit to significant sites from linear model
  linear_results <- read.csv(paste0('GTEx/GTEx',tissue,'_RNAseq_linear_results.csv'))
  sig_genes <- linear_results[linear_results$fdr < 0.1,]$molecule
  if(length(sig_genes) < 1000) {
    sig_genes <- linear_results[order(linear_results$fdr),]$molecule[1:1000]
  }
  print(length(sig_genes))
  tissue_data <- tissue_data[rownames(tissue_data) %in% sig_genes,]
  dim(tissue_data)
  print(identical(colnames(tissue_data), rownames(tissue_metadata)))
  timepoints <- tissue_metadata$AGE
  print(nrow(tissue_metadata))
  
  filename <- paste0("GTEx/GTEx",tissue,"_RNAseq_GTEx_lmer_results.csv")
  datafilename <- paste0("GTEx/GTEx",tissue,"_RNAseq_data.RData")
  metadatafilename <-paste0('GTEX/GTEx_',tissue,'RNAseq_metadata.csv')
  #save(tissue_data,file=datafilename)
  # write.csv(tissue_metadata,metadatafile)
  cohort <- as.factor(tissue_metadata$COHORT)
  sex <- as.factor(tissue_metadata$SEX)
  race <- as.factor(tissue_metadata$RACE)
  ethnicity <- as.factor(tissue_metadata$ETHNCTY)
  bmi <- as.numeric(tissue_metadata$BMI)
  tissue_metadata$MHSMKSTS[tissue_metadata$MHSMKSTS == ""] <- 'No'
  smoker <- as.factor(tissue_metadata$MHSMKSTS)
  random_effects <- data.frame(cohort,sex,race,ethnicity,bmi,smoker)
  
  print("modeling..")
  results <- GTEx_lmer_function(tissue_data,timepoints,random_effects,filename)
  head(results)
  dim(results)
  nSubTissues <- length(unique(tissue_metadata$Note))
  model <- 'LMER'
  fileData[i,] <- c(tissue, nrow(tissue_metadata),filename,nSubTissues,datafilename,metadatafilename,model)
}
write.csv(fileData,"GTEx/GTEx_TissueFile_GTEx_lmer_metadata.csv")

##########

## Linear Model for GTEx MultiTissue RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
subTissues <- metadata$Note
table(subTissues)

# Break up by tissue
tissues <- unique(subTissues)
nTissues <- length(tissues)
nTissues
fileData <- matrix(ncol=8,nrow=nTissues)
colnames(fileData) <- c("Tissue","nSamples","filename","nSubTissues","datafilename","metadatafilename","model","model_parameters")
for(i in 1:nTissues) {
  tissue <- tissues[i]
  print(i)
  print(tissue)
  tissue_metadata <- metadata[metadata$Note == tissue,]
  tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
  print(identical(colnames(tissue_data), rownames(tissue_metadata)))
  timepoints <- tissue_metadata$AGE
  print(nrow(tissue_metadata))
  if(tissue == "Cells - Leukemia cell line (CML)") { # nrow(tissue_metadata) < 50 |
    next
  }
  print("modeling..")
  filename <- paste0("GTEx/linear/GTEx",tissue,"_RNAseq_linear_results.csv")
  datafilename <- paste0("GTEx/data/GTEx",tissue,"_RNAseq_data.RData")
  #save(tissue_data,file=datafilename)
  metadatafile <- paste0('GTEX/metadata/GTEx_',tissue,'RNAseq_metadata.csv')
  write.csv(tissue_metadata,metadatafile)
  # results <- linear_function(tissue_data,timepoints,filename)
  # head(results)
  # dim(results)
  nSubTissues <- length(unique(tissue_metadata$Note))
  model <- 'linearRegression'
  model_parameters <- ""
  fileData[i,] <- c(tissue, nrow(tissue_metadata),filename,nSubTissues,datafilename,metadatafile,model,model_parameters)
}
write.csv(fileData,"GTEx/GTEx_TissueFile_linear_metadata.csv")

##########

## MinMax Model for GTEx MultiTissue RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
subTissues <- metadata$Note
table(subTissues)

# Break up by tissue
tissues <- unique(subTissues)
nTissues <- length(tissues)
nTissues
fileData <- matrix(ncol=8,nrow=nTissues)
colnames(fileData) <- c("Tissue","nSamples","filename","nSubTissues","datafilename","metadatafilename","model","model_parameters")
for(i in 1:nTissues) {
  tissue <- tissues[i]
  
  print(tissue)
  tissue_metadata <- metadata[metadata$Note == tissue,]
  tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
  print(identical(colnames(tissue_data), rownames(tissue_metadata)))
  timepoints <- tissue_metadata$AGE
  print(nrow(tissue_metadata))
  if(tissue == "Cells - Leukemia cell line (CML)" | nrow(tissue_metadata) < 50 ) { # 
    next
  }
  print("modeling..")
  filename <- paste0("GTEx/minmax/GTEx",tissue,"_RNAseq_minmax_results.csv")
  datafilename <- paste0("GTEx/data/GTEx",tissue,"_RNAseq_data.RData")
  #save(tissue_data,file=datafilename)
  metadatafile <- paste0('GTEX/metadata/GTEx_',tissue,'RNAseq_metadata.csv')
  write.csv(tissue_metadata,metadatafile)
  # results <- MinMax_function(tissue_data,timepoints,filename)
  # head(results)
  # dim(results)
  nSubTissues <- length(unique(tissue_metadata$Note))
  model <- 'confidenceIntervals'
  model_parameters <- ""
  fileData[i,] <- c(tissue, nrow(tissue_metadata),filename,nSubTissues,datafilename,metadatafile,model,model_parameters)
}
write.csv(fileData,"GTEx/GTEx_TissueFile_minmax_metadata.csv")

##########

## Permutation Model for GTEx MultiTissue RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
subTissues <- metadata$Note
table(subTissues)

# Break up by tissue
tissues <- unique(subTissues)
nTissues <- length(tissues)
nTissues
fileData <- matrix(ncol=8,nrow=nTissues)
colnames(fileData) <- c("Tissue","nSamples","filename","nSubTissues","datafilename","metadatafilename","model","model_parameters")
for(i in 1:nTissues) {
  print(i)
  tissue <- tissues[i]
  print(tissue)
  if(tissue == "Cells - Leukemia cell line (CML)") {
    next
  }
  tissue_metadata <- metadata[metadata$Note == tissue,]
  tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
  # limit to significant sites from linear model
  linear_results <- read.csv(paste0('GTEx/linear/GTEx',tissue,'_RNAseq_linear_results.csv'))
  sig_genes <- linear_results[linear_results$fdr < 0.1,]$molecule
  if(length(sig_genes) < 1000) {
    sig_genes <- linear_results[order(linear_results$fdr),]$molecule[1:1000]
  }
  print(length(sig_genes))
  tissue_data <- tissue_data[rownames(tissue_data) %in% sig_genes,]
  dim(tissue_data)
  print(identical(colnames(tissue_data), rownames(tissue_metadata)))
  timepoints <- tissue_metadata$AGE
  print(nrow(tissue_metadata))

  filename <- paste0("GTEx/permutation/GTEx",tissue,"_RNAseq_permutation_results.csv")
  datafilename <- paste0("GTEx/data/GTEx",tissue,"_RNAseq_data.RData")
  metadatafile <- paste0('GTEX/metadata/GTEx_',tissue,'RNAseq_metadata.csv')
  nSubTissues <- length(unique(tissue_metadata$Note))
  model <- 'permutation'
  # nBatchSampling <- round(nrow(tissue_metadata)*0.1,0)
  model_parameters <- 'pc0.01_ns100_nb1000'
  if(nrow(tissue_metadata) <= 50) {
    next
  }
  # save(tissue_data,file=datafilename)
  # if(nBatchSampling < 5) {
  #   next
  # }
  fileData[i,] <- c(tissue, nrow(tissue_metadata),filename,nSubTissues,datafilename,metadatafile,model,model_parameters)
  
  print("modeling..")
  
  results <- Permutation_function(tissue_data,timepoints,filename,nBootstrap,nBatchSampling,percentile)
  head(results)
  dim(results)
  
}
write.csv(fileData,"GTEx/GTEx_TissueFile_permutation_metadata.csv")
##########

## Permutation Model for GTEx Skin RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
subTissues <- metadata$Note
tissue <- 'Skin - Sun Exposed (Lower leg)'
tissue_metadata <- metadata[metadata$Note == tissue,]
tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
# limit to significant sites from linear model
linear_results <- read.csv('GTEx/GTExSkin - Sun Exposed (Lower leg)_RNAseq_linear_results.csv')
sig_genes <- linear_results[linear_results$fdr < 0.05,]$molecule
length(sig_genes)
tissue_data <- tissue_data[rownames(tissue_data) %in% sig_genes,]
dim(tissue_data)
print(identical(colnames(tissue_data), rownames(tissue_metadata)))
timepoints <- tissue_metadata$AGE
print(nrow(tissue_metadata))
filename <- paste0("GTEx/GTExSkin - Sun Exposed (Lower leg)_RNAseq_permutation_pc0.01_ns100_nb1000_results.csv")
print("modeling..")
results <- Permutation_function(tissue_data,timepoints,filename,nBooststrap,nBatchSampling,percentile)
head(results)
dim(results)

##########

## GTEx_LMER Model for GTEx Skin RNAseq
##########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
subTissues <- metadata$Note
tissue <- 'Skin - Sun Exposed (Lower leg)'
tissue_metadata <- metadata[metadata$Note == tissue,]
tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
# limit to significant sites from linear model
linear_results <- read.csv(paste0('GTEx/linear/GTEx',tissue,'_RNAseq_linear_results.csv'))
sig_genes <- linear_results[linear_results$fdr < 0.1,]$molecule
print(length(sig_genes))
tissue_data <- tissue_data[rownames(tissue_data) %in% sig_genes,]
dim(tissue_data)
print(identical(colnames(tissue_data), rownames(tissue_metadata)))
timepoints <- tissue_metadata$AGE
cohort <- as.factor(tissue_metadata$COHORT)
sex <- as.factor(tissue_metadata$SEX)
race <- as.factor(tissue_metadata$RACE)
ethnicity <- as.factor(tissue_metadata$ETHNCTY)
bmi <- as.numeric(tissue_metadata$BMI)
tissue_metadata$MHSMKSTS[tissue_metadata$MHSMKSTS == ""]<- 'No'
smoker <- as.factor(tissue_metadata$MHSMKSTS )
random_effects <- data.frame(cohort,sex,race,ethnicity,bmi,smoker)

filename <- paste0("GTEx/GTEx_lmer/GTExSkin - Sun Exposed (Lower leg)_RNAseq_GTEx_lmer_results.csv")
print("modeling..")
results <- GTEx_lmer_function(tissue_data,timepoints,random_effects,filename)
head(results)
dim(results)

##########

##### Residual Variance of GTEx Fibroblasts ########
load("GTEx/GTEx_RNAseq_data.RData")
GTEx_RNAseq_data[1:5,1:5]
max(GTEx_RNAseq_data)
min(GTEx_RNAseq_data)
dim(GTEx_RNAseq_data)
metadata <- read.csv("GTEx/GTEx_RNAseq_metadata.csv",row.names=1)
metadata[1:5,1:3]
identical(colnames(GTEx_RNAseq_data),rownames(metadata))
tissue <- metadata$Tissue
table(tissue)
tissue <- 'Skin'
tissue_metadata <- metadata[metadata$Tissue == tissue,]
tissue_data <-GTEx_RNAseq_data[,colnames(GTEx_RNAseq_data) %in% rownames(tissue_metadata)]
print(identical(colnames(tissue_data), rownames(tissue_metadata)))
timepoints <- tissue_metadata$AGE
skin_results <- read.csv('GTEx/GTExSkin_RNAseq_lmer_results.csv')
skin_results <- skin_results[order(skin_results$fdr,decreasing=T),]
skin_results$molecule[1:10]
gene <- skin_results$molecule[1]
top_data <- tissue_data[rownames(tissue_data) == gene,]
model_data <- data.frame(timepoints,top_data)
model <- lm(top_data ~ timepoints,data=model_data)
intercept <- coef(model)[1]
intercept
slope <- coef(model)[2]
slope
sum <- summary(model)
r2 <- sum$r.squared
r2
pValue <- sum$coefficients[2,4]
pValue
residuals <- residuals(model)
residuals
model_summary <- c(gene, 'linearRegression',slope, intercept, pValue, r2)
names(model_summary) <- c('gene','model_type','slope','intercept','pValue','r2')
model_results <- data.frame(colnames(tissue_data),model_data, residuals)
colnames(model_results) <- c('SamName','Age','Expression','Residual')
hist(residuals,breaks=50)
p <- ggplot(data=model_data,aes(x=timepoints,y=top_data)) +
  geom_point(size=1,alpha=0.5) +
  #stat_summary(fun.data=mean_cl_normal) +
  scale_y_continuous(name=gene) +
  geom_smooth(method='lm', formula= y~x)
setwd(baseDir)
write.csv(model_summary,paste0('GTEx/GTEx_Skin_',gene,'linear_model_summary.csv'))
write.csv(model_results,paste0('GTEx/GTEx_Skin_',gene,'linear_model_residuals.csv'))
pdf(paste('GTEx/GTEx_Skin_',gene,'_plots.pdf'))
p
hist(residuals,breaks=50)
dev.off()
################
          
######## PART 2 - RESCALING ########

 # Variables: Molecule, Molecular_Annotation, Dataset, nSamples, Species, System, StudyName, Assay, Platform, TissueType , Disease, Treatment, FittedModel, slope, intercept, pValue, fdr, r2, 
Species <- c('Human','Mouse')
Systems <- c('Primary_Cell_Culture','Biopsy','PostMortem','Immortalized_Cell_Culture')
Studies <- c('Cellular_Lifespan', 'TwinsUk', 'iPOP','NIA','GTEx','Tabula_Muris')
Tissues <- c('Fibroblasts','Fat','Subcutaneous_Fat', 'Brown_Fat', 
             'Skin','Muscle','Bone','Marrow','Whole_Blood','White_Blood_Cells',
             'Brain','Heart','Liver','Pancreas','Spleen','Lung','Kidney',)
Assays <- c('Transcriptomics','DNAmethylation','Proteomics','Metabolomics')
Platforms <- c('RNAseq','RNAarray','EPICarray','450karray')
Treatments <- c('Control','Contact_Inhibition','Hypoxia','Oligomycin')
Diseases <- c('Healthy','SURF1_Mutation')
FittedModels <- c('linearRegression','LMER','Permuataion','MinMax')


###### Aggregate into Single Matrix ####
### Load in Dataset Slopes ###
setwd(baseDir)
dir()
cols <- c("Dataset","comparison","nSamples", "Species", "System", "StudyName", "TissueType" , "Assay", "Platform",  "Treatment", "Disease" , "model","filename","metadata")
datasets <- list()
addMetdata <- function(ds, specs) {
  specframe <- data.frame(matrix(ncol=length(cols),nrow=nrow(ds)))
  for(i in 1:length(specs)) { 
    specframe[,i] <- rep(specs[i], length(nrow(ds)))
  }
  colnames(specframe) <- cols
  ds <- data.frame(ds, specframe)
  datasets <- list.append(datasets, ds)
  return(datasets)
}

##### DNAm #####
main_DNAm_comparison <- 'TwinsUK_Skin_DNAm_linearRegression'
#main_DNAm_comparison <- 'TwinsUK_Skin_DNAm_permutation'
ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_DNAm','all',38,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_metadata.csv')
datasets <- addMetdata(ds,specs)

# ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_linear_linear_results.csv")
# specs <- c('Control_Fibroblasts_r1_linear_DNAm',main_DNAm_comparison,23,'Human',
#            'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
#            'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
#            'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_linear_betas.RData',
#            'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_linear_metadata.csv')
# datasets <- addMetdata(ds,specs)
# 
# ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_senescent_linear_results.csv")
# specs <- c('Control_Fibroblasts_r1_senescent_DNAm',main_DNAm_comparison,15,'Human',
#            'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
#            'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
#            'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescent_betas.RData',
#            'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescent_metadata.csv')
# datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_earlyPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_earlyPhase_DNAm',main_DNAm_comparison,12,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_midPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_midPhase_DNAm',main_DNAm_comparison,11,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_senescentPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_senescentPhase_DNAm',main_DNAm_comparison,12,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_treatments_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_treatments_DNAm',main_DNAm_comparison,15,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r2_linear_results.csv")
specs <- c('Control_Fibroblasts_r2_DNAm',main_DNAm_comparison,23,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r3_linear_results.csv")
specs <- c('Control_Fibroblasts_r3_DNAm',main_DNAm_comparison,18,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/contact_inhibition_fibroblasts_linear_results.csv")
specs <- c('Contact_Inhibition_Fibroblasts_DNAm','all',21,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/contact_inhibition_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Contact_Inhibition_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,9,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/contact_inhibition_latePhase_fibroblasts_linear_results.csv")
specs <- c('Contact_Inhibition_latePhase_Fibroblasts_DNAm',main_DNAm_comparison,12,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_latePhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/DEX_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('DEX_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,15,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','100nM dexamethasone','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/SURF1_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('SURF1_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,14,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','control','SURF1_mutation','linearRegression',
           'Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/modulators_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Modulators_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,14,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','mitoNUITs','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/oligo_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Oligomycin_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,13,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','mitoNUITs','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_linear_results.csv")
specs <- c('TwinsUK_Skin_DNAm','all',322,'Human','Biopsy','TwinsUk',
           'Skin','DNAmethylation','450karray','Control','Healthy','linearRegression',
           'TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData',
           'TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_permutation_results.csv")
specs <- c('TwinsUK_Skin_DNAm','all',322,'Human','Biopsy','TwinsUk',
           'Skin','DNAmethylation','450karray','Control','Healthy','permutation',
           'TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData',
           'TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv')
datasets <- addMetdata(ds,specs)


# ds <- read.csv("TwinsUK/DNAm_Fat/twinsUK_Fat_DNAm_linear_results.csv")
# specs <- c('TwinsUK_Fat_DNAm','all',648,'Human','Biopsy','TwinsUk',
#            'Fat','DNAmethylation','450karray','Control','Healthy','linearRegression',
#            'TwinsUK/DNAm_Fat/twinsUK_FAT_DNAm_betas.RData',
#            'TwinsUK/DNAm_Fat/RNAarray_age_values.csv')
# datasets <- addMetdata(ds,specs)

# ds <- read.csv("Horvath/UCLA_fibroblasts_DNAm_results.csv")
# specs <- c('UCLA_Fibroblasts_DNAm','all',147,'Human','Biopsy','UCLA','Fibroblasts','DNAmethylation','450karray','Control','Healthy','biweightMidcorrelation')
# datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/HEK293_DNAm_linear_results.csv")
specs <- c('HEK293_DNAm','all',12,'Human','Immortalized_Cell_Culture','Cellular_Lifespan',
           'Embryonic_Kidney_Cells','DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/HEK293_betas.RData',
           'Cellular_Lifespan/DNAm/HEK293_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/Hypoxia_fibroblasts_linear_results.csv")
specs <- c('Hypoxia_Fibroblasts_DNAm',main_DNAm_comparison,26,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Hypoxia_3%_O2','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Hypoxia_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/deoxyglucose_fibroblasts_linear_results.csv")
specs <- c('Deoxyglucose_Fibroblasts_DNAm',main_DNAm_comparison,11,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','2-Deoxyglucose','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)
#########

##### RNAseq #####
main_RNAseq_comparison <- 'GTEx_Skin - Sun Exposed (Lower leg)_RNAseq_linearRegression'
#main_RNAseq_comparison <- 'GTEx_Skin - Sun Exposed (Lower leg)_RNAseq_permutation'
ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_RNAseq','all',29,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

# ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_linear_results.csv")
# specs <- c('Control_Fibroblasts_r1_linear_RNAseq',main_RNAseq_comparison,15,'Human','Primary_Cell_Culture',
#            'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
#            'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_genes.RData',
#            'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_linear_RNAseq_metadata.csv')
# datasets <- addMetdata(ds,specs)
# 
# ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_linear_results.csv")
# specs <- c('Control_Fibroblasts_r1_senescent_RNAseq',main_RNAseq_comparison,14,'Human','Primary_Cell_Culture',
#            'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
#            'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_genes.RData',
#            'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescent_RNAseq_metadata.csv')
# datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_earlyPhase_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_midPhase_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_senescentPhase_RNAseq',main_RNAseq_comparison,11,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_treatments_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r2_RNAseq',main_RNAseq_comparison,16,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r3_RNAseq',main_RNAseq_comparison,19,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Contact_Inhibition_Fibroblasts_RNAseq','all',21,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/contact_inhibition_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/contact_inhibition_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Contact_Inhibition_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Contact_Inhibition_latePhase_Fibroblasts_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/contact_inhibition_latePhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Hypoxia_Fibroblasts_RNAseq',main_RNAseq_comparison,15,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Hypoxia_3%_O2','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/hypoxia_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Hypoxia_Contact_Inhibition_Fibroblasts_RNAseq',main_RNAseq_comparison,15,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq',
           'Hypoxia_3%_O2_Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/hypoxia_contact_inhibition_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_linear_results.csv")
specs <- c('DEX_Fibroblasts_RNAseq',main_RNAseq_comparison,30,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Dexamethasone','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/dex_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('DEX_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Dexamethasone','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_linear_results.csv")
specs <- c('SURF1_Fibroblasts_RNAseq',main_RNAseq_comparison,20,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','SURF1_Mutation','linearRegression',
           'Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/surf1_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('SURF1_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','SURF1_Mutation','linearRegression',
           'Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Oligomycin_Fibroblasts_RNAseq',main_RNAseq_comparison,20,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Oligomycin','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/oligo_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/oligo_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Oligomycin_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,10,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','UK5099+BPTES+Etomoxir','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Mito_Modulators_Fibroblasts_RNAseq',main_RNAseq_comparison,23,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','UK5099+BPTES+Etomoxir','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/modulators_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Mito_Modulators_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,11,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','UK5099+BPTES+Etomoxir','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_linear_results.csv")
specs <- c('2-Deoxyglucose_Fibroblasts_RNAseq',main_RNAseq_comparison,24,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','2-Deoxyglucose','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/deoxyglucose_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('2-Deoxyglucose_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','2-Deoxyglucose','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Betahydroxybutyrate_Fibroblasts_RNAseq',main_RNAseq_comparison,19,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Betahydroxybutyrate','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/betahydroxybutyrate_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Galactose_Fibroblasts_RNAseq',main_RNAseq_comparison,18,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Galactose','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/galactose_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

# ds <- read.csv("iPOP/iPOP_RNAseq_lmer_results.csv")
# specs <- c('iPOP_Whole_Blood_RNAseq','Control_Fibroblasts_r1_RNAseq_linearRegression',881,'Human','Biopsy',
#            'iPOP','Whole_Blood','Transcriptomics','RNAseq','Control','Healthy','LMER',
#            '')
# datasets <- addMetdata(ds,specs)

# ds <- read.csv("NIA/NIA_RNAseq_linear_results.csv")
# specs <- c('NIA_Muscle_RNAseq','Control_Fibroblasts_r1_RNAseq_linearRegression',53,'Human','Biopsy',
#            'NIA','Muscle','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
#            '')
# datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_linear_results.csv")
specs <- c('TwinsUK_Skin_RNAarray','all',705,'Human','Biopsy',
           'TwinsUK','Skin','Transcriptomics','RNAarray','Control','Healthy','linearRegression',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_permutation_results.csv")
specs <- c('TwinsUK_Skin_RNAarray','all',705,'Human','Biopsy',
           'TwinsUK','Skin','Transcriptomics','RNAarray','Control','Healthy','permutation',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/RNAarray/twinsUK_skin_RNAarray_minmax_results.csv")
specs <- c('TwinsUK_Skin_RNAarray','all',705,'Human','Biopsy',
           'TwinsUK','Skin','Transcriptomics','RNAarray','Control','Healthy','confidenceIntervals',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_genes.RData',
           'TwinsUK/RNAarray/twinsUK_skin_RNAarray_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/RNAarray/twinsUK_fat_RNAarray_linear_results.csv")
specs <- c('TwinsUK_Fat_RNAarray','all',825,'Human','Biopsy',
           'TwinsUK','Fat','Transcriptomics','RNAarray','Control','Healthy','linearRegression',
           'TwinsUK/RNAarray/twinsUK_fat_RNAarray_genes.RData',
           'TwinsUK/RNAarray/twinsUK_fat_RNAarray_metadata.csv')
datasets <- addMetdata(ds,specs)

# ds <- read.csv("TwinsUK/RNAarray/twinsUK_adipose_RNAarray_lmer_results.csv")
# specs <- c('TwinsUK_Adipose_RNAarray','all',825,'Human','Biopsy','TwinsUK','Adipose','Transcriptomics','RNAarray','Control','Healthy','LMER')
# datasets <- addMetdata(ds,specs)

# ds <- read.csv("TwinsUK/RNAarray/twinsUK_lymphocytes_RNAarray_lmer_results.csv")
# specs <- c('TwinsUK_Lymphocytes_RNAarray','Control_Fibroblasts_r1_RNAseq',92,'Human','Biopsy',
#            'TwinsUK','Lymphocytes','Transcriptomics','RNAarray','Control','Healthy','LMER',
#            '')
# datasets <- addMetdata(ds,specs)

# ds <- read.csv("TwinsUK/RNAarray/twinsUK_LCLs_RNAarray_lmer_results.csv")
# specs <- c('TwinsUK_LCLs_RNAarray','Control_Fibroblasts_r1_RNAseq',825,'Human','Biopsy',
#            'TwinsUK','LCLs','Transcriptomics','RNAarray','Control','Healthy','LMER',
#            '')
# datasets <- addMetdata(ds,specs)

ds <- read.csv("Tabula_Muris/Tabula_Muris_RNAseq_lmer_results.csv")
specs <- c('Tabula_Muris_MultiTissue_RNAseq','all',892,'Mouse','Biopsy',
           'TwinsUK','MultiTissue','Transcriptomics','RNAseq','Control','Healthy','LMER',
           'Tabula_Muris/Tabula_Muris_RNAseq_genes.RData',
           'Tabula_Muris/Tabula_Muris_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

GTExfiles <- read.csv("GTEx/GTEx_TissueFile_linear_metadata.csv",header=T)
for(i in 1:nrow(GTExfiles)) {
  tissue <- GTExfiles$Tissue[i]
  if(is.na(tissue)) {next}
  filename <- GTExfiles$filename[i]
  datafilename <- GTExfiles$datafilename[i]
  metadatafilename <- GTExfiles$metadatafilename[i]
  modelType <- GTExfiles$model[i]
  nSamples <- GTExfiles$nSamples[i]
  comparison <- 'Control_Fibroblasts_r1_RNAseq_linearRegression'
  if(tissue == 'Skin - Sun Exposed (Lower leg)' | tissue == 'Cells - Transformed fibroblasts') {
    comparison <- "all"
  }
  ds <- read.csv(filename)
  dataset_name <- paste0('GTEx_',tissue,'_RNAseq') #  
  specs <- c(dataset_name,comparison,nSamples,'Human','PostMortem','GTEx',tissue,'Transcriptomics','RNAseq','Control','Healthy',modelType,datafilename,metadatafilename)
  datasets <- addMetdata(ds,specs)
}

GTExfiles <- read.csv("GTEx/GTEx_TissueFile_minmax_metadata.csv",header=T)
for(i in 1:nrow(GTExfiles)) {
  tissue <- GTExfiles$Tissue[i]
  if(is.na(tissue)) {next}
  filename <- GTExfiles$filename[i]
  datafilename <- GTExfiles$datafilename[i]
  metadatafilename <- GTExfiles$metadatafilename[i]
  modelType <- GTExfiles$model[i]
  nSamples <- GTExfiles$nSamples[i]
  comparison <- 'Control_Fibroblasts_r1_RNAseq_linearRegression'
  # if(tissue == 'Skin - Sun Exposed (Lower leg)' | tissue == 'Cells - Transformed fibroblasts') {
  #   comparison <- "all"
  # }
  ds <- read.csv(filename)
  dataset_name <- paste0('GTEx_',tissue,'_RNAseq')
  specs <- c(dataset_name,comparison,nSamples,'Human','PostMortem','GTEx',tissue,'Transcriptomics','RNAseq','Control','Healthy',modelType,datafilename,metadatafilename)
  datasets <- addMetdata(ds,specs)
}

GTExfiles <- read.csv("GTEx/GTEx_TissueFile_permutation_metadata.csv",header=T)
for(i in 1:nrow(GTExfiles)) {
  tissue <- GTExfiles$Tissue[i]
  if(is.na(tissue)) {next}
  filename <- GTExfiles$filename[i]
  datafilename <- GTExfiles$datafilename[i]
  metadatafilename <- GTExfiles$metadatafilename[i]
  modelType <- GTExfiles$model[i]
  nSamples <- GTExfiles$nSamples[i]
  comparison <- 'Control_Fibroblasts_r1_RNAseq_linearRegression'
  if(tissue == 'Skin - Sun Exposed (Lower leg)' | tissue == 'Cells - Transformed fibroblasts') {
    comparison <- "all"
  }
  ds <- read.csv(filename)
  dataset_name <- paste0('GTEx_',tissue,'_RNAseq')
  specs <- c(dataset_name,comparison,nSamples,'Human','PostMortem','GTEx',tissue,'Transcriptomics','RNAseq','Control','Healthy',modelType,datafilename,metadatafilename)
  datasets <- addMetdata(ds,specs)
}


ds <- read.csv("GTEx/GTEx_lmer/GTExSkin - Sun Exposed (Lower leg)_RNAseq_GTEx_lmer_results.csv")
specs <- c('GTEx_Skin_RNAseq','all',520,'Human','Biopsy',
           'GTEx','Skin - Sun Exposed (Lower leg)','Transcriptomics','RNAseq','Control','Healthy','LMER',
           'GTEx/data/GTExSkin - Sun Exposed (Lower leg)_RNAseq_data.RData', 
           'GTEx/metadata/GTEx_Skin - Sun Exposed (Lower leg)RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)
#############

##### DNAm - treatments #####
main_DNAm_comparison <- 'TwinsUK_Skin_DNAm_linearRegression'
#main_DNAm_comparison <- 'TwinsUK_Skin_DNAm_permutation'
ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_DNAm','all',38,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r2_linear_results.csv")
specs <- c('Control_Fibroblasts_r2_DNAm',main_DNAm_comparison,23,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r2_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r3_linear_results.csv")
specs <- c('Control_Fibroblasts_r3_DNAm',main_DNAm_comparison,18,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r3_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_earlyPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_earlyPhase_DNAm',main_DNAm_comparison,12,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_earlyPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_midPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_midPhase_DNAm',main_DNAm_comparison,11,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_midPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_senescentPhase_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_senescentPhase_DNAm',main_DNAm_comparison,12,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_senescentPhase_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/control_fibroblasts_r1_treatments_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_treatments_DNAm',main_DNAm_comparison,15,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Control','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_betas.RData',
           'Cellular_Lifespan/DNAm/Control_Fibroblasts_r1_treatments_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/contact_inhibition_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Contact_Inhibition_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,9,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Contact_Inhibition_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/DEX_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('DEX_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,15,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','100nM dexamethasone','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/DEX_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/SURF1_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('SURF1_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,14,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','control','SURF1_mutation','linearRegression',
           'Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/SURF1_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/modulators_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Modulators_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,14,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','mitoNUITs','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Modulators_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <-  read.csv("Cellular_Lifespan/DNAm/oligo_earlyPhase_fibroblasts_linear_results.csv")
specs <- c('Oligomycin_earlyPhase_Fibroblasts_DNAm',main_DNAm_comparison,13,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','mitoNUITs','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Oligo_earlyPhase_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_linear_results.csv")
specs <- c('TwinsUK_Skin_DNAm','all',322,'Human','Biopsy','TwinsUk',
           'Skin','DNAmethylation','450karray','Control','Healthy','linearRegression',
           'TwinsUK/DNAm_Skin/twinsUK_Skin_DNAm_betas.RData',
           'TwinsUK/DNAm_Skin/DNAm_GEO_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/DNAm/deoxyglucose_fibroblasts_linear_results.csv")
specs <- c('Deoxyglucose_Fibroblasts_DNAm',main_DNAm_comparison,11,'Human',
           'Primary_Cell_Culture','Cellular_Lifespan','Fibroblasts',
           'DNAmethylation','EPICarray','2-Deoxyglucose','Healthy','linearRegression',
           'Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_betas.RData',
           'Cellular_Lifespan/DNAm/Deoxyglucose_Fibroblasts_metadata.csv')
datasets <- addMetdata(ds,specs)
#########

##### RNAseq -treatments #####
main_RNAseq_comparison <- 'GTEx_Skin - Sun Exposed (Lower leg)_RNAseq_linearRegression'
#main_RNAseq_comparison <- 'GTEx_Skin - Sun Exposed (Lower leg)_RNAseq_permutation'
ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_RNAseq','all',29,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_earlyPhase_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_earlyPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r2_RNAseq',main_RNAseq_comparison,16,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r2_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r3_RNAseq',main_RNAseq_comparison,19,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r3_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_midPhase_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_midPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_senescentPhase_RNAseq',main_RNAseq_comparison,11,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_senescentPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_linear_results.csv")
specs <- c('Control_Fibroblasts_r1_treatments_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/control_fibroblasts_r1_treatments_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Contact_Inhibition_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture','Cellular_Lifespan',
           'Fibroblasts','Transcriptomics','RNAseq','Contact_Inhibition','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/contact_inhibition_earlyPhase_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('DEX_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Dexamethasone','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/dex_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('SURF1_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,12,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','Control','SURF1_Mutation','linearRegression',
           'Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/surf1_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/oligo_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Oligomycin_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,10,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','UK5099+BPTES+Etomoxir','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/oligomycin_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('Mito_Modulators_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,11,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','UK5099+BPTES+Etomoxir','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/modulators_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

ds <- read.csv("Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_linear_results.csv")
specs <- c('2-Deoxyglucose_earlyPhase_Fibroblasts_RNAseq',main_RNAseq_comparison,9,'Human','Primary_Cell_Culture',
           'Cellular_Lifespan','Fibroblasts','Transcriptomics','RNAseq','2-Deoxyglucose','Healthy','linearRegression',
           'Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_genes.RData',
           'Cellular_Lifespan/RNAseq/deoxyglucose_earlyPhase_fibroblasts_RNAseq_metadata.csv')
datasets <- addMetdata(ds,specs)

GTExfiles <- read.csv("GTEx/GTEx_TissueFile_linear_metadata.csv",header=T)
for(i in 1:nrow(GTExfiles)) {
  tissue <- GTExfiles$Tissue[i]
  if(is.na(tissue)) {next}
  filename <- GTExfiles$filename[i]
  datafilename <- GTExfiles$datafilename[i]
  metadatafilename <- GTExfiles$metadatafilename[i]
  modelType <- GTExfiles$model[i]
  nSamples <- GTExfiles$nSamples[i]
  comparison <- 'Control_Fibroblasts_r1_RNAseq_linearRegression'
  if(tissue == 'Skin - Sun Exposed (Lower leg)' | tissue == 'Cells - Transformed fibroblasts') {
    comparison <- "all"
  }
  ds <- read.csv(filename)
  dataset_name <- paste0('GTEx_',tissue,'_RNAseq') #  
  specs <- c(dataset_name,comparison,nSamples,'Human','PostMortem','GTEx',tissue,'Transcriptomics','RNAseq','Control','Healthy',modelType,datafilename,metadatafilename)
  datasets <- addMetdata(ds,specs)
}
#############

###### Compare Datasets #######
setwd(baseDir)
cutoff <- 0.1
minMolecules <- 2000
plotGraphs <- F
GAMrestrict <- F
SIGrestrict <- F
rescaling_results_filename <- "rescaling_results/summary_6_indNonSig.csv"

# Main Graphing Function
rescale_function <- function(dataset1, dataset2, filename, graphPlots, GAM, SIG) {
  colnames(dataset2) <- colnames(dataset1)
  dataset1 <- dataset1[order(dataset1$fdr),]
  dataset2 <- dataset2[order(dataset2$fdr),]
  dataset1_sig <- dataset1
  dataset2_sig <- dataset2
  
  if(SIG == T) {
    # Significant Molecules
    dataset1_sig <- dataset1[dataset1$fdr <= cutoff,]
    if(nrow(dataset1_sig) < minMolecules) {
      dataset1_sig <- dataset1[1:minMolecules,]
    }
    dataset2_sig <- dataset2[dataset2$fdr <= cutoff,]
    if(nrow(dataset2_sig) < minMolecules) {
      dataset2_sig <- dataset2[1:minMolecules,]
    }
  }
  
  if(GAM == T) {
    GAM_control_results <- matrix()
    if(dataset1$Assay[1] == "DNAmethylation") {
      setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Robert_Analysis/")
      GAM_control_results <- read.csv("DNAm_GAM_sig_sites_3controls.csv")
    }
    else if(dataset1$Assay[1] == "Transcriptomics") {
      genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
      setwd(genewizDir)
      GAM_control_results <- read.csv("RNAseq_GAM_sig_sites_3controls.csv")
    }
    setwd(baseDir)
    ### filter to GAM sig sites
    dataset1_sig <- dataset1_sig[dataset1_sig$molecule %in% GAM_control_results$X,]
    dataset2_sig <- dataset2_sig[dataset2_sig$molecule %in% GAM_control_results$X,]
  }
  
  

  # Median of abs slopes
  dataset1_median <- median(abs(dataset1_sig$slope), na.rm = TRUE)
  dataset2_median <- median(abs(dataset2_sig$slope), na.rm = TRUE)
  
  # if dataset2 is dominant to dataset1 flip datasets
  if(dataset2_median > dataset1_median) {
    holddata <- dataset1
    dataset1 <- dataset2
    dataset2 <- holddata
    
    holddata <- dataset1_median
    dataset1_median <- dataset2_median
    dataset2_median <- holddata
    
    holddata <- dataset1_sig
    dataset1_sig <- dataset2_sig
    dataset2_sig <- holddata
  }
  
  all_sig_data <- rbind(dataset1_sig,dataset2_sig)
  #print(head(all_sig_data))
  
  study1 <- paste0(dataset1$Dataset[1],'_',dataset1$model[1])
  print(study1)
  study2 <- paste0(dataset2$Dataset[1],'_',dataset2$model[1])
  print(study2)
  
  molecule = "genes"
  assay = "Transcriptomics"
  if(dataset1_sig$Assay[1] == 'DNAmethylation') {
    molecule = "CpGs"
    assay="DNAmethylation"
  }
  # Systems <- c('Primary_Cell_Culture','Biopsy','PostMortem','Immortalized_Cell_Culture')
  system1 <- "InVitro"
  if(dataset1$System[1] == "Biopsy" | dataset1$System[1] == "PostMortem") {
    system1 <- "InVivo"
  }
  system2 <- "InVitro"
  if(dataset2$System[1] == "Biopsy" | dataset2$System[1] == "PostMortem") {
    system2 <- "InVivo"
  }
  comparison_type <- paste0(system1,"-",system2)
  print(comparison_type)
  
  ### Distribution of slopes
  print("slopes...")
  medians <- c(dataset2_median,dataset1_median)
  fc_pos <- max(medians) / 2
  median_fc <- max(medians) / min(medians)
  Slope_distribution <- ggplot(data=all_sig_data, aes(x=abs(slope), color=Dataset,fill=Dataset)) +
    geom_histogram(aes(y=..density..),alpha=1, position ="identity", binwidth = 0.01, size=0.3) + # 
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = dataset1_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset1_median,2),y=0,x=dataset1_median), vjust=-1,col='red',size=5)+
    geom_vline(xintercept = dataset2_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset2_median,2),y=0,x=dataset2_median), vjust=-1,col='red',size=5)+
    geom_text(aes(label=paste0("Slope Fold Change: ",signif(median_fc,3)),y=0,x=fc_pos), vjust=-10,col='red',size=5)+
    # 
    scale_fill_manual(values=c("green", "purple")) +
    scale_color_manual(values=c("green", "purple")) +
    scale_x_log10(name = "Slope (log10, fdr<=0.05)") + # limits=c(-20,0)
    scale_y_continuous(name = "Density") +
    ggtitle(label= paste0("Distribution of slopes (all significant ",molecule,")")) +
    coord_cartesian(clip = 'off') + 
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  #print(Slope_distribution)
  
  ### Breakdown of sig molecules
  
  # join results
  dataset1_sig_overlap <- dataset1_sig[dataset1_sig$molecule %in% dataset2_sig$molecule,]
  dataset2_sig_overlap <- dataset2_sig[dataset2_sig$molecule %in% dataset1_sig_overlap$molecule,]
  dataset1_sig_overlap <- dataset1_sig_overlap[match(dataset2_sig_overlap$molecule,dataset1_sig_overlap$molecule),]
  print(identical(dataset1_sig_overlap$molecule, dataset2_sig_overlap$molecule))
  overlap_data <- data.frame(dataset1_sig_overlap, dataset2_sig_overlap)
  print(dim(overlap_data))
  rescaling_factor <- overlap_data$slope / overlap_data$slope.1
  rescaling_factor_abs <- abs(overlap_data$slope) / abs(overlap_data$slope.1)
  combined_pvalue <- -log(overlap_data$pValue) + -log(overlap_data$pValue.1)
  combined_delta_pvalue <-combined_pvalue - (-log(overlap_data$pValue) - -log(overlap_data$pValue.1))
  overlap_data <- data.frame(overlap_data,rescaling_factor,rescaling_factor_abs, combined_pvalue, combined_delta_pvalue)
  dim(overlap_data)
  
  
  ### Distribution of overlapping slopes
  print("slopes2...")
  dataset1_overlap_median <- median(abs(dataset1_sig_overlap$slope), na.rm = TRUE)
  dataset2_overlap_median <- median(abs(dataset2_sig_overlap$slope), na.rm = TRUE)
  all_overlap_data <- rbind(dataset1_sig_overlap, dataset2_sig_overlap)
  medians <- c(dataset2_overlap_median,dataset1_overlap_median)
  fc_pos <- max(medians) / 2
  median_overlap_fc <- max(medians) / min(medians)
  Slope_distribution_overlap <- ggplot(data=all_overlap_data, aes(x=abs(slope), color=Dataset,fill=Dataset)) +
    geom_histogram(aes(y=..density..),alpha=1, position ="identity", binwidth = 0.01, size=0.3) + # 
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = dataset1_overlap_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset1_overlap_median,2),y=0,x=dataset1_overlap_median), vjust=-1,col='red',size=5)+
    geom_vline(xintercept = dataset2_overlap_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset2_overlap_median,2),y=0,x=dataset2_overlap_median), vjust=-1,col='red',size=5)+
    geom_text(aes(label=paste0("Slope Fold Change: ",signif(median_overlap_fc,3)),y=0,x=fc_pos), vjust=-10,col='red',size=5)+
    scale_fill_manual(values=c("green", "purple")) +
    scale_color_manual(values=c("green", "purple")) +
    scale_x_log10(name = "Slope (log10, fdr<=0.05)") + # limits=c(-20,0)
    scale_y_continuous(name = "Density") +
    ggtitle(label= "Distribution of overlapping slopes (both directions)") +
    coord_cartesian(clip = 'off') + 
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  
  study <- factor(c(study1,study2, "Total",study1,study2, "Total",study1,study2, "Total"), levels = c(study1,study2, "Total"))
  total_genes <- c(nrow(dataset1_sig), nrow(dataset2_sig), nrow(overlap_data)) 
  combined_upregulated <- nrow(overlap_data[overlap_data$slope > 0 & overlap_data$slope.1 > 0,])
  combined_upregulated
  combined_downregulated <- nrow(overlap_data[overlap_data$slope < 0 & overlap_data$slope.1 < 0,])
  combined_downregulated
  combined_updown <- nrow(overlap_data) - combined_upregulated - combined_downregulated
  combined_updown
  up_genes <- c(nrow(dataset1_sig[dataset1_sig$slope > 0,]), nrow(dataset2_sig[dataset2_sig$slope>0,]), combined_upregulated)
  down_genes <- c(nrow(dataset1_sig[dataset1_sig$slope < 0,]), nrow(dataset2_sig[dataset2_sig$slope<0,]), combined_downregulated)
  mixed_genes <- c(0,0, combined_updown)
  genes <- c(up_genes, down_genes, mixed_genes)
  label = factor(c(rep("upregulated",3), rep("downregulated",3), rep("mixed",3)), levels=c("upregulated","mixed","downregulated"))
  if(dataset1_sig$Assay[1] == 'DNAmethylation') {
    label = factor(c(rep("Hypermethylated",3), rep("Hypomethylated",3), rep("mixed",3)),levels=c("Hypermethylated","mixed","Hypomethylated"))
  }
  totaldata <- data.frame(study, genes, label)
  head(totaldata)
  print("barplot...")
  sig_genes_barplot <- ggplot(totaldata, aes(x = study, y=genes, fill = label)) +
    geom_bar(stat="identity", colour="black") +
    scale_fill_manual(values =c('darkred','darkgray','darkblue')) + 
    geom_text(aes(label=genes), size=3.5,color="white",position = position_stack(vjust = 0.5)) + 
    ggtitle(label= paste("Breakdown of significant ",molecule)) +
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  print("hist...")
  pos_mode <- median(overlap_data[overlap_data$rescaling_factor > 0,]$rescaling_factor,na.rm=T)
  print(pos_mode)
  neg_mode <-median(overlap_data[overlap_data$rescaling_factor <  0,]$rescaling_factor,na.rm=T)
  print(neg_mode)
  maxx <- pos_mode * 5
  minx <- -maxx
  Rescaling_factor_hist <- ggplot(data=overlap_data, aes(x=rescaling_factor)) +
    geom_histogram(aes(y = ..density..), alpha=1, position ="identity", bins =200, color = "grey30", fill = "gray", size=0.3) +#bins <- 500 # binwidth = 0.01,
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = pos_mode, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(pos_mode,2),y=0,x=pos_mode),
              vjust=-5,col='red',size=4)+
    geom_vline(xintercept = neg_mode, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(neg_mode,2),y=0,x=neg_mode),
              vjust=-5,col='red',size=4)+
    #scale_fill_manual(values=c("green", "purple")) +
    #scale_color_manual(values=c("green", "purple")) +
    scale_x_continuous(name = "Rescaling Factor",limits=c(minx,maxx)) + # limits=c(-20,0)
    scale_y_continuous(name = "Density") + # limits=c(-2000,2000)
    ggtitle(label= "Distribution of Rescaling Factor (fdr<0.05)") +
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  Rescaling_factor_hist
  
  print("hist2...")
  median_abs_rf <- median(overlap_data$rescaling_factor_abs, na.rm=T)
  maxx <- median_abs_rf *5
  minx <- 0
  Rescaling_factor_hist_abs <- ggplot(data=overlap_data, aes(x=rescaling_factor_abs)) +
    geom_histogram(aes(y = ..density..), alpha=1, position ="identity", bins =200, color = "grey30", fill = "gray", size=0.3) +#bins <- 500 # binwidth = 0.01,
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = median_abs_rf, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(median_abs_rf,2),y=0,x=median_abs_rf),
              vjust=-7,col='red',size=5)+
    #scale_fill_manual(values=c("green", "purple")) +
    #scale_color_manual(values=c("green", "purple")) +
    scale_x_continuous(name = "Rescaling Factor (abs)",limits=c(minx,maxx)) + # limits=c(-20,0)
    scale_y_continuous(name = "Density") + # limits=c(-2000,2000)
    ggtitle(label= "Distribution of Rescaling Factor (fdr<0.05), absolute value of slopes") +
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  Rescaling_factor_hist_abs
  
  print("hist3...")
  overlap_data_sd <- overlap_data[overlap_data$rescaling_factor >0,]
  median_rf <- median(overlap_data_sd$rescaling_factor, na.rm=T)
  log2_median_rf <- log2(median_rf)
  maxx <- log2_median_rf *2
  minx <- 0
  Rescaling_factor_hist_sd <- ggplot(data=overlap_data_sd, aes(x=log2(rescaling_factor))) +
    geom_histogram(aes(y = ..density..), alpha=1, position ="identity", bins =200, color = "grey30", fill = "gray", size=0.3) +#bins <- 500 # binwidth = 0.01,
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = log2_median_rf, color = "red", linetype = "dashed") +
    geom_text(aes(label=paste0(signif(log2_median_rf,2),'\n',signif(median_rf,2)),y=0,x=log2_median_rf),
              vjust=-7,col='red',size=5)+
    #scale_fill_manual(values=c("green", "purple")) +
    #scale_color_manual(values=c("green", "purple")) +
    scale_x_continuous(name = "Rescaling Factor (log2)",limits=c(minx,maxx)) + # limits=c(-20,0)
    scale_y_continuous(name = "Density") + # limits=c(-2000,2000)
    ggtitle(label= "Distribution of Rescaling Factor (fdr<0.05), only same direction slopes") +
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  Rescaling_factor_hist_sd
  
  ### Distribution of overlapping slopes in the same direction
  print("slopes3...")
  dataset1_sig_overlap_sd <- dataset1_sig_overlap[dataset1_sig_overlap$molecule %in% overlap_data_sd$molecule,]
  dataset2_sig_overlap_sd <- dataset2_sig_overlap[dataset2_sig_overlap$molecule %in% overlap_data_sd$molecule,]
  
  dataset1_overlap_sd_median <- median(abs(dataset1_sig_overlap_sd$slope), na.rm = TRUE)
  dataset2_overlap_sd_median <- median(abs(dataset2_sig_overlap_sd$slope), na.rm = TRUE)
  all_overlap_sd_data <- rbind(dataset1_sig_overlap_sd, dataset2_sig_overlap_sd)
  medians <- c(dataset2_overlap_sd_median,dataset1_overlap_sd_median)
  fc_pos <- max(medians) / 2
  median_overlap_sd_fc <- max(medians) / min(medians)
  print(paste0("# of overlap same direction molecules: ",nrow(all_overlap_sd_data)))
  Slope_distribution_overlap_sd <- ggplot(data=all_overlap_sd_data, aes(x=abs(slope), color=Dataset,fill=Dataset)) +
    geom_histogram(aes(y=..density..),alpha=1, position ="identity", binwidth = 0.01, size=0.3) + # 
    geom_density(alpha = .2, fill = "antiquewhite3") + 
    geom_vline(xintercept = dataset1_overlap_sd_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset1_overlap_sd_median,2),y=0,x=dataset1_overlap_sd_median), vjust=-1,col='red',size=5)+
    geom_vline(xintercept = dataset2_overlap_sd_median, color = "red", linetype = "dashed") +
    geom_text(aes(label=signif(dataset2_overlap_sd_median,2),y=0,x=dataset2_overlap_sd_median), vjust=-1,col='red',size=5)+
    geom_text(aes(label=paste0("Slope Fold Change: ",signif(median_overlap_sd_fc,3)),y=0,x=fc_pos), vjust=-10,col='red',size=5)+
    scale_fill_manual(values=c("green", "purple")) +
    scale_color_manual(values=c("green", "purple")) +
    scale_x_log10(name = "Slope (log10, fdr<=0.05)") + # limits=c(-20,0)
    scale_y_continuous(name = "Density") +
    ggtitle(label= "Distribution of overlapping slopes (same direction)") +
    coord_cartesian(clip = 'off') + 
    theme_classic(base_size = 16) +
    theme(axis.line = element_line(colour = 'black', size = 0.5),
          plot.title = element_text(size = 12, face = "bold"))
  
  print("Scatterplot...")
  Slopes_scatter <- ggplot(overlap_data_sd, aes(x=slope, y = slope.1)) +
    geom_point(size = 1, alpha = 0.5) +
    theme_classic(base_size = 20) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    scale_x_continuous(name=paste0(study1," slope (delta/year)")) +
    scale_y_continuous(name=paste0(study2, " slope (delta/year)")) +
    ggtitle(label= "Scatterplot of slopes, same direction") +
    theme(axis.line = element_line(colour = 'black', size = 0.2),
          plot.title = element_text(size = 12, face = "bold"))
 
  overlap_data_sd <- overlap_data_sd[order(overlap_data_sd$combined_delta_pvalue,decreasing=T),]
  data_filename <- paste0("rescaling_results/data/",study1,"_",study2,"overlap_sd.csv")
  write.csv(overlap_data_sd, data_filename)
  print('timecourses...')
  if(graphPlots == T) {
    # Plot top 10 sites
    load.Rdata(file=dataset1$filename[1], 'betas1')
    print(betas1[1:3,1:3])
    print(dim(betas1))
    metadata1 <- read.csv(dataset1$metadata[1])
    load.Rdata(file=dataset2$filename[1], 'betas2')
    print(betas2[1:3,1:3])
    print(dim(betas2))
    metadata2 <- read.csv(dataset2$metadata[1])
    timepoints1 <- ""
    if('Timepoints_years' %in% colnames(metadata1)) {
      timepoints1 <- metadata1$Timepoints_years
    }
    else if ('age' %in% colnames(metadata1)) {
      timepoints1 <- metadata1$age
    }
    else if ('Age' %in% colnames(metadata1)) {
      timepoints1 <- metadata1$Age
    }
    else if ('AGE' %in% colnames(metadata1)) {
      timepoints1 <- metadata1$AGE
    }
    else if ('age.at.biopsy.ch1' %in% colnames(metadata1)) {
      timepoints1 <- metadata1$age.at.biopsy.ch1
    }
    
    
    timepoints2 <- ""
    if('Timepoints_years' %in% colnames(metadata2)) {
      timepoints2 <- metadata2$Timepoints_years
    }
    else if ('age' %in% colnames(metadata2)) {
      timepoints2 <- metadata2$age
    }
    else if ('Age' %in% colnames(metadata2)) {
      timepoints2 <- metadata2$Age
    }
    else if ('AGE' %in% colnames(metadata2)) {
      timepoints2 <- metadata2$AGE
    }
    else if ('age.at.biopsy.ch1' %in% colnames(metadata2)) {
      timepoints2 <- metadata2$age.at.biopsy.ch1
    }
    timepoints <- c(timepoints1,timepoints2)
    
    subjects1 <- rep("",length(timepoints1))
    if('Cell_Line' %in% colnames(metadata1)) {
      subjects1 <- metadata1$Cell_Line
    }
    subjects2 <- rep("",length(timepoints2))
    if('Cell_Line' %in% colnames(metadata2)) {
      subjects2 <- metadata2$Cell_Line
    }
    subjects = c(subjects1,subjects2)
    plot_top_sites_ctrl <- function(overlap_data_sd,betas1,betas2,metadata1,metadata2,timepoints,subjects,n_sites) {
      model1 <- overlap_data_sd$model[1]
      model2 <- overlap_data_sd$model.1[1]
      top_sites <- overlap_data_sd$molecule[1:n_sites]
      for(i in 1:length(top_sites)) {
        Molecule <- top_sites[i]
        #Molecule <- "NID1"
        molecule_data1 <- betas1[rownames(betas1) == Molecule,]
        molecule_data2 <- betas2[rownames(betas2) == Molecule,]
        molecule_data <- c(molecule_data1,molecule_data2)
        studyID<-c(rep(study1,length(molecule_data1)),rep(study2,length(molecule_data2)))
        Site_Data <- data.frame(molecule_data, studyID, subjects, timepoints)
        yintercept1 <- signif(overlap_data_sd$intercept[i],2)
        yintercept2 <- signif(overlap_data_sd$intercept.1[i],2)
        slope1 <- signif(overlap_data_sd$slope[i],2)
        slope2 <- signif(overlap_data_sd$slope.1[i],2)
        pValue1 <- signif(overlap_data_sd$pValue[i],2)
        pValue2 <- signif(overlap_data_sd$pValue.1[i],2)
        r21 <- signif(overlap_data_sd$r2[i],2)
        r22 <- signif(overlap_data_sd$r2.1[i],2)
        anno <- paste0('\n ',study1,'\n slope: ',slope1,' intercept: ',yintercept1,' pValue: ',pValue1,' r2: ',r21,
                       '\n',study2,'\n slope: ',slope2,' intercept: ',yintercept2,' pValue: ',pValue2,' r2: ',r22)
        p <- ggplot(Site_Data, aes(x=timepoints, y = molecule_data, color = studyID)) +
          geom_point(size = 2 , stroke = 1, aes(fill = studyID, shape = subjects)) +
          geom_abline(intercept = yintercept1, slope = slope1, color="green", size=1)+
          geom_abline(intercept = yintercept2, slope = slope2, color="purple", size=1)+
          #annotate(geom="text", x=0, y=1, label=anno, color="black") +
          scale_color_manual(values = c("green", "purple")) +
          scale_fill_manual(values =  c("green","purple")) +
          #scale_shape_manual(values=c(21,22,23,21)) +
          scale_y_continuous(name = Molecule) +
          scale_x_continuous(name = "Age (years)") + 
          theme_classic() +
          ggtitle(label= paste(Molecule, " #",i,anno,sep="")) +
          coord_cartesian(clip = "off") +
          theme(legend.position='bottom') +
          facet_grid(cols = vars(studyID),scales = "free")
        # if(Log == T) {
        #   p = p + scale_x_log10(name = "Log10 Age (years)")
        # }
        # else {
        #   p = p + scale_x_continuous(name = "Age (years)")
        # }
        print(p)
      }
    }
    # plot_top_sites_ctrl(overlap_data_sd,betas1,betas2,metadata1,metadata2,timepoints,cell_lines, 5, F)
  }

  
  print('tables...')
  col1 <- as.vector(t(dataset1[1,8:ncol(dataset1)]))
  #print(col1)
  col2 <- as.vector(t(dataset2[1,8:ncol(dataset2)]))
  #print(col2)
  Table <- data.frame(col1,col2)
  colnames(Table) <- c(study1,study2)
  #print(Table)
  rownames(Table) <- colnames(dataset1)[8:ncol(dataset1)]
  nSignificant_Molecules <- c(nrow(dataset1_sig), nrow(dataset2_sig))
  Table <- rbind(Table, nSignificant_Molecules)
  rownames(Table)[15] <- paste0('# of sig. ',molecule)
  # add median slopes
  summary_slopes <- data.frame(c(dataset1_median, dataset1_overlap_median, dataset1_overlap_sd_median),
                               c(dataset2_median, dataset2_overlap_median, dataset2_overlap_sd_median))
  colnames(summary_slopes) <- c(study1,study2)
  Table <- rbind(Table, summary_slopes)
  rownames(Table)[16:18] <- c(paste0('median slope all sig. ',molecule),paste0('median slope of shared ',molecule,' (both directions)'),paste0('median slope of shared ',molecule,' (same direction)'))
  Table2 <- c(paste(study1,study2,sep="\n"),assay,comparison_type,round(median_fc,2),
              round(median_overlap_fc,2),round(median_overlap_sd_fc,2),
              nrow(overlap_data),nrow(overlap_data_sd),
              combined_upregulated,combined_downregulated,combined_updown)
  names(Table2) <- c('Groups','Assay','Comparison_Type',
                     paste0('Rescaling Factor of all sig. ',molecule),
                     paste0('Rescaling Factor of shared ',molecule,' (both directions)'),
                     paste0('Rescaling Factor of shared ',molecule,' (same direction)'),
                     paste0('# of shared sig. ',molecule),
                     paste0('# of shared sig. ',molecule,' (same direction)'),
                     paste0('# of upregulated ',molecule),
                     paste0('# of downregualted ',molecule),
                     paste0('# of mixed direction ',molecule))
  if(graphPlots == T) {
    print("Graphing...")
    ### Save Graphs
    setwd(baseDir)
    pdf(filename)
    grid.table(Table, rows = rownames(Table),theme=ttheme_default(base_size = 6))
    grid.newpage()
    grid.table(Table2, rows = names(Table2),theme=ttheme_default(base_size = 6))
    #grid.newpage()
    print(sig_genes_barplot)
    print(Slope_distribution)
    print(Slope_distribution_overlap)
    print(Slope_distribution_overlap_sd)
    print(Rescaling_factor_hist)
    print(Rescaling_factor_hist_abs)
    print(Rescaling_factor_hist_sd)
    print(Slopes_scatter)
    
    # hist(filtered_fbr_betas, main="Fibroblasts Controls DNAm Histogram", xlab = "beta", breaks = 100)
    # hist(fat_betas, main="TwinsUK DNAm Histogram", xlab = "beta", breaks = 100)
    plot_top_sites_ctrl(overlap_data_sd,betas1,betas2,metadata1,metadata2,timepoints,subjects, 5)
    dev.off()
    rm(betas1,betas2)
  }

  
  # Organize results into single row
  r1 <- Table2
  r1[1] <- paste(study1,study2,sep=" vs ")
  names(r1) <- c('group','assay','comparison','rf1','rf2','rf3','nShared','nSD','nUP','nDown','mMixed')
  r2 <- Table[,1]
  names(r2) <- paste0("ds1_",rownames(Table))
  r3 <- Table[,2]
  names(r3) <- paste0("ds2_",rownames(Table))
  results <- c(r1,r2,r3)
  print("")
  return(results)
}

nDatasets <- length(datasets)
nDatasets
rescale_results <- vector()
count <- 0
for(i in 1:(nDatasets-1)) {
  dataset1 <- datasets[[i]]
  for(j in (i+1):nDatasets) {
    count <- count +1
    #print(count)
    dataset2 <- datasets[[j]]
    if(dataset1$Assay[1] == dataset2$Assay[1]) {
      if(dataset1$comparison[1] != 'all') {
        if(dataset1$comparison[1] != paste0(dataset2$Dataset[1],"_",dataset2$model[1])) {
          # print('break1')
          next
        }
      }
      if(dataset2$comparison[1] != 'all') {
        if(dataset2$comparison[1] != paste0(dataset1$Dataset[1],"_",dataset1$model[1])) {
          #print('break2')
          next
        }
      }
      dataset1_sig <- dataset1[dataset1$fdr <= cutoff,]
      dataset2_sig <- dataset2[dataset2$fdr <= cutoff,]
      # if(nrow(dataset1_sig) < minMolecules | nrow(dataset2_sig) < minMolecules) {
      #   print('break3')
      #   next
      # }
      filename <- paste("rescaling_results/graphs/rescaled_graphs",dataset1$Dataset[1],dataset2$Dataset[1],sep="_")
      filename <- paste0(filename,".pdf")
      print(count)
      rescale_results <- rbind(rescale_results, rescale_function(dataset1,dataset2,filename,plotGraphs,GAMrestrict,SIGrestrict))
    }
  }
}
write.csv(rescale_results,rescaling_results_filename)


#### Summary Plots #####
# Graph summary rescaling results bar plot
rescale_results <- as.data.frame(rescale_results)
rescaling_factor <- c(rescale_results$rf1,rescale_results$rf2,rescale_results$rf3)
groups <- gsub(" v. ", "\n", rep(rescale_results$group,3))
assay <- rep(rescale_results$assay,3)
comparison_type <- rep(rescale_results$comparison,3)
rescaling_type <- c(rep("rf1",nrow(rescale_results)), 
                    rep("rf2",nrow(rescale_results)), 
                    rep("rf3",nrow(rescale_results)))
rescale_results_df <- data.frame(groups,assay,comparison_type,rescaling_factor,rescaling_type)

rescaling_barplot <- ggplot(data=rescale_results_df,aes(x=groups, y=as.numeric(rescaling_factor), fill=rescaling_type)) +
  geom_bar(stat = "identity",position=position_dodge(),size=.3,alpha=1,color='black') + 
  geom_text(aes(label=round(as.numeric(rescaling_factor),1)), position=position_dodge(width=0.9),vjust=-0.5,size=5) +
  #coord_flip() +
  scale_y_log10(name="Resacling Factor (log10)")+
  ggtitle(label= "Comparison of Rescaling Factors")+
  scale_fill_manual(values=c("chartreuse3","darkgoldenrod2","chocolate2"),
                    labels = c("Rescaling Factor (all sig. molecules)", "Rescaling Factor (shared both directions)", "Rescaling Factor (shared same direction)")) +
  coord_cartesian(clip = 'off') + 
  theme_classic(base_size = 30) +
  facet_wrap(~assay+comparison_type,scales = "free_x",ncol=3) +
  theme(axis.line = element_line(colour = 'black', size = 0.2),
        axis.text.x = element_text(angle = 45,hjust=1.001),
        axis.title.y = element_text(margin = margin(t = 0, r = 50, b = 0, l = 0)),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position="bottom",legend.title = element_blank())
#rescaling_barplot

pdf('rescaling_results/summary_barplot_4.pdf', width=100,height=25)
rescaling_barplot
dev.off()

# Correlate sample size with rescaling factor
deltaSampleSize <-  abs(as.numeric(rescale_results$ds2_nSamples) - as.numeric(rescale_results$ds1_nSamples))
groups <- unique(rescale_results$comparison)
sampleSizePlot <- ggplot(rescale_results, aes(x=as.numeric(deltaSampleSize),y=as.numeric(rf3), color=comparison)) +
  geom_point(size=2,alpha=0.5) +
  geom_smooth(method='lm', formula= y~x,se =F) +
  stat_poly_eq(formula = y~x, aes(group = comparison, label =paste(groups, stat(eq.label), stat(rr.label), sep = "~~")),
                       rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "top", show.legend = TRUE) +
  #geom_text(label=rescale_results[,1],hjust=-.05,size=1) +
  scale_y_log10(name="Log10 Rescaling Factor (same direction)") + # limits=c(0,1000)
  scale_x_continuous(name="delta sample size") +
  ggtitle(label= "Correlation of delta sample size and rescaling factor (same direction)")+
  coord_cartesian(clip = 'off') + 
  theme_classic(base_size = 12) +
  theme(axis.line = element_line(colour = 'black', size = 0.2),
        plot.title = element_text(size = 10, face = "bold"),
        legend.position="right",legend.title = element_blank())
#sampleSizePlot

pdf('rescaling_results/summary_sampleSize.pdf', width=6,height=6)
sampleSizePlot
dev.off()
###############
  
