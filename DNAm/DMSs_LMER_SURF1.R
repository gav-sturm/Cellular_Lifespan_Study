#####################################
##
## Title: DNA Methylation Linear Mixed Effects Modeling to obtain differntially methylated CpGs (DMSs)
## Author: Gabriel Sturm 
## Date: 2020-08-01
## see lme4 documentation: https://cran.r-project.org/web/packages/lme4/lme4.pdf
##
##
##


library(minfi)                                         
library(limma)                                         
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) # annotation for Illumina's EPIC methylation arrays
library(ggplot2)
library(mgcv)
library(lme4)


baseDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/"
setwd(baseDir)


baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")

#load lifespan study data
baseDir3 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp"
setwd(baseDir3)
dir()
LS_Data <- read.csv("Lifespan_Study_Data.csv")


setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/")
#write.csv(IlluminaAnnot.EPIC, "Illumina_EPIC_Annotation.csv")
IlluminaAnnot.EPIC <- read.csv("Illumina_EPIC_Annotation.csv") # 

# Filter betas to sites in annotation
IlluminaAnnot.EPIC <-  IlluminaAnnot.EPIC[order(match(IlluminaAnnot.EPIC$Name,rownames(betas))),]  
betas <- betas[rownames(betas) %in% IlluminaAnnot.EPIC$Name,]
IlluminaAnnot.EPIC <- IlluminaAnnot.EPIC[IlluminaAnnot.EPIC$Name %in% rownames(betas),]



# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14","hFB6","hFB7","hFB8")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatments == "Control",]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2,]
filtered_data <-  filtered_data[filtered_data$basename != "",]

# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown > 15,]
filtered_data <- filtered_data[filtered_data$Days_Grown < 75,]

# Split data by control vs SURF1 samples
Control_data <- filtered_data[filtered_data$Clinical_Condition == "Normal",]
SURF1_data <- filtered_data[filtered_data$Clinical_Condition == "SURF1_Mutation",]

# Retreive betas for control and beta samples
filtered_betas <- betas[,colnames(betas) %in% as.vector(filtered_data$basename)]
Control_betas <- betas[,colnames(betas) %in% as.vector(Control_data$basename)]
SURF1_betas <- betas[,colnames(betas) %in% as.vector(SURF1_data$basename)]

# Order betas in same order as Lifespan Data
filtered_betas <- filtered_betas[,match(filtered_data$basename, colnames(filtered_betas))]
Control_betas <- Control_betas[,match(Control_data$basename, colnames(Control_betas))]
SURF1_betas <- SURF1_betas[,match(SURF1_data$basename, colnames(SURF1_betas))]

Control_timepoints <- Control_data$Days_Grown
Control_timepoints
length(Control_timepoints)
SURF1_timepoints <- SURF1_data$Days_Grown
SURF1_timepoints
length(SURF1_timepoints)

# Source: https://rcompanion.org/handbook/I_09.html
# Assumptions - The ANOVA test makes the following assumptions about the data:
# 1. Independence of the observations. Each subject should belong to only one group. There is no relationship between the observations in each group. Having repeated measures for the same participants is not allowed.
# 2. No significant outliers in any cell of the design
# 3. Normality. the data for each design cell should be approximately normally distributed.
# 4. Homogeneity of variances. The variance of the outcome variable should be equal in every cell of the design.
LMER_test <- function(filtered_betas, filtered_data) {
  P_values <- matrix(ncol=1, nrow=nrow(filtered_betas))
  for(i in 1:nrow(filtered_betas)) {
    site_betas <- filtered_betas[i,]
    # Test Sites
    #site_betas <- filtered_betas[rownames(filtered_betas) == "cg21959717",]
    #site_betas <- filtered_betas[rownames(filtered_betas) == "cg16248866",]
    CpG <- data.frame(filtered_data, site_betas)
    
   
    # Mixed Effects Model
    # memodel = lmer(site_betas ~ Clincal_Condition*Days_Grown + (1+Days_Grown|Cell_Line), data = CpG,  REML = FALSE)
    # summary(memodel)
    # P_values[i,] <- summary(memodel)$coefficients[,5]
    # 
    # SURF1_Data <- subset_samples[subset_samples$Treatments == "Control",]
    # SURF1_betas <- subset_betas[,colnames(subset_betas) %in% SURF1_Data$basename]
    
    MEmodel = lmer(site_betas ~ Clinical_Condition*Days_Grown + (1|Days_Grown),
                   data = CpG,  REML = F)
    summary <- summary(MEmodel)
    
    null_model <- lmer(site_betas ~ Days_Grown + (1+Days_Grown|Cell_Line),
                       data = CpG,  REML = F)
    
    anova <- anova(null_model, MEmodel)
    P_values[i,1] <- anova$`Pr(>Chisq)`[2]
    # print(ggplot(CpG, aes(x=Days_Grown, y=site_betas,
    #                              shape=Cell_Line,color=Clinical_Condition)) + geom_point() + geom_smooth(method=lm, se=F))
    message(i)
  }
  # Correct for multipule Comparisons
  bonfi_P_values <- as.matrix(p.adjust(P_values[,1], method = "bonferroni", n = length(P_values)))
  fdr_P_values <- as.matrix(p.adjust(P_values[,1], method = "fdr", n = length(P_values)))
  
  #Organize matrix
  results <- data.frame(P_values, bonfi_P_values, fdr_P_values)
  rownames(results) <- rownames(filtered_betas)
  colnames(results) <- c("P_values", "Bonfi_P_values", "FDR_P_values")
  return(results)
}


results <- LMER_test(filtered_betas, filtered_data)
# setwd(baseDir)
# write.csv(results, "LMER_SURF1_50_days.csv")

### Plot the top sites of results
#order sites by category
results <- results[order(results$P_values),]


#age_results <- Group_Results_Annotation[order(Group_Results_Annotation$methyl_difference, decreasing = TRUE),]
# choose number of sites to plot
plot_top_sites <- function(n_sites, results, linear) {
  top_sites <- rownames(results)[1:n_sites]
  for(i in 1:length(top_sites)) {
    CpG <- top_sites[i]
    #CpG <- "cg16867657"
    Site_betas <- c(Control_betas[rownames(Control_betas) == CpG,],SURF1_betas[rownames(SURF1_betas) == CpG,]) 
    Site_Data <- cbind(filtered_data,Site_betas)
    
    p <- ggplot(Site_Data, aes(x=Days_Grown, y = Site_betas)) +
      geom_point(size = 2 , stroke = 0.5, aes(shape = Cell_Line, fill = Cell_Line)) +
      #geom_line(aes(y = Control_GAM)) +
      scale_color_manual(values = c("dimgray", "blue")) +
      scale_fill_manual(values =  c("dimgray","dimgray","dimgray","blue","blue","blue")) +
      scale_shape_manual(values=c(21,22,23,21,22,23)) +
      scale_y_continuous(name = CpG, breaks = seq(0, 1.0, by = 0.2), limits= c(0,1.0)) +
      scale_x_continuous(name =  "Time Grown (Days)", breaks = seq(0, 300, by = 50), limits = c(0,300)) +
      theme_classic() +
      coord_cartesian(clip = "off") +
      theme(legend.position='none')
    
    if(linear == TRUE) {
      p = p + geom_smooth(aes(group = Cell_Line, col = Clinical_Condition),method = lm, se = FALSE)
      cell_lines <- Site_Data$Cell_Line
      # p = p + stat_poly_eq(formula = y~x, aes(group = cell_lines, label =  paste(cell_lines, stat(eq.label), stat(rr.label), sep = "~~")),
      #                                             rr.digits = 2, coef.digits = 2,parse = TRUE,  label.x = "left", label.y = "bottom", show.legend = TRUE)
    }
    else {
      control_gam <- gam(Control_betas[rownames(betas) == CpG,]~s(Control_timepoints,bs="cr",k=6))
      SURF1_gam <- gam(SURF1_betas[rownames(betas) == CpG,]~s(SURF1_timepoints,bs="cr",k=6)) 
      #   use predict(gam) & summary(gam) to plot all info
      c_pv <- round(summary(control_gam)$s.pv,2)
      c_r2 <- round(summary(control_gam)$r.sq,2)
      c_dof <- round(summary(control_gam)$edf,2)
      
      s_pv <- round(summary(SURF1_gam)$s.pv,2)
      s_r2 <- round(summary(SURF1_gam)$r.sq,2)
      s_dof <- round(summary(SURF1_gam)$edf,2)
      
      annotation <- paste("GAM RESULTS \n Control: \n Pvalue = ", c_pv, "\n R2 = ", c_r2,"\n DoF = ", c_dof, "\n \n SURF1: \n Pvalue = ", s_pv, "\n R2 = ", s_r2,"\n DoF = ", s_dof, sep = "")
      
      p = p + geom_smooth(aes(group = Cell_Line, col = Clinical_Condition),method = gam, formula = y ~ s(x, k=6, bs = "cr"), se = FALSE)
      p = p + annotate("text", label = annotation, x = 0.5, y = 0.5, size = 3, hjust = 0, vjust = 1)
    }
    print(p)
  }
}
plot_top_sites(10, results, FALSE)

# Function for caluating the methylation differnce between two lists of betas
methyl_difference <- function(group1_betas, group2_betas) {
  # Split Sites by hyper and hypo methylation SURF1 vs control
  Methyl_diff <- matrix(ncol=2, nrow=nrow(group1_betas))
  for(i in 1:nrow(group1_betas)) {
    diff <- mean(group2_betas[i,]) - mean(group1_betas[i,])
    Methyl_diff[i,2] <- diff
    if(diff > 0) {
      Methyl_diff[i,1] <- "Hypermethylated"
    }
    else {
      Methyl_diff[i,1] <- "Hypomethylated"
    }
  }
  rownames(Methyl_diff) <- rownames(group1_betas)
  colnames(Methyl_diff) <- c("methyl_direction", "methyl_difference")
  return(Methyl_diff)
}


# Run methyl differnce comparing control vs SURF1
methyl_diff <- methyl_difference(Control_betas, SURF1_betas) # 1min
head(methyl_diff)
methyl_diff <-  methyl_diff[order(match(rownames(methyl_diff),rownames(results))),]  

results <- data.frame(results, methyl_diff)
head(results)

# Attach results to annotation
IlluminaAnnot.EPIC <-  IlluminaAnnot.EPIC[order(match(IlluminaAnnot.EPIC$Name,rownames(results))),]  
results_in_anno <- results[rownames(results) %in% IlluminaAnnot.EPIC$Name,]
IlluminaAnnot.EPIC <- IlluminaAnnot.EPIC[IlluminaAnnot.EPIC$Name %in% rownames(results_in_anno),]
Results_Annotation <- cbind(results_in_anno, IlluminaAnnot.EPIC)
Results_Annotation <- Results_Annotation[order(Results_Annotation$Bonfi_P_value),]
head(Results_Annotation)

cutoff <- 0.01
sig_results_anno <- Results_Annotation[Results_Annotation$FDR_P_values < cutoff,]
nrow(sig_results_anno) # 5,028

setwd(baseDir)
write.csv(sig_results_anno, "LMER_SURF1_50_days_sig_anno.csv")

