#####################################
##
## Title: DNA Methylation CombP to obtain differntially methylated regions (DMRs)
## Author: Gabriel Sturm 
## contact info: gs2934@cumc.columbia.edu
## Date: 2020-10-01
## see ENmix documentation: https://rdrr.io/bioc/ENmix/f/inst/doc/ENmix.pdf
##
##
##

#BiocManager::install("ENmix")
library(ENmix)
library(doParallel)
library(ggplot2)

baseDir2 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/Preprocessing/data/"
setwd(baseDir2)
dir()
load("Combat_Betas.RData")


setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/")
IlluminaAnnot.EPIC <- read.csv("Illumina_EPIC_Annotation.csv") # 

# Filter betas to sites in annotation
IlluminaAnnot.EPIC <-  IlluminaAnnot.EPIC[order(match(IlluminaAnnot.EPIC$Name,rownames(betas))),]  
betas <- betas[rownames(betas) %in% IlluminaAnnot.EPIC$Name,]
IlluminaAnnot.EPIC <- IlluminaAnnot.EPIC[IlluminaAnnot.EPIC$Name %in% rownames(betas),]

tss15 <- IlluminaAnnot.EPIC[IlluminaAnnot.EPIC$corrected_gene_group == 'TSS1500',]

#load lifespan study data
baseDir3 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp"
setwd(baseDir3)
dir()
LS_Data <- read.csv("Lifespan_Study_Data.csv")


# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14","hFB6","hFB7","hFB8")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatments == "Control",]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2,]
filtered_data <-  filtered_data[filtered_data$basename != "",]
filtered_data$basename
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown > 15,]
filtered_data <- filtered_data[filtered_data$Days_Grown < 75,]

# Split data by control vs SURF1 samples
Control_data <- filtered_data[filtered_data$Clinical_Condition == "Normal",]
SURF1_data <- filtered_data[filtered_data$Clinical_Condition == "SURF1_Mutation",]

# Timepoints
Control_timepoints <- Control_data$Days_Grown
Control_timepoints
length(Control_timepoints)
SURF1_timepoints <- SURF1_data$Days_Grown
SURF1_timepoints
length(SURF1_timepoints)

# Retreive betas for control and beta samples
filtered_betas <- betas[,colnames(betas) %in% as.vector(filtered_data$basename)]
Control_betas <- betas[,colnames(betas) %in% as.vector(Control_data$basename)]
SURF1_betas <- betas[,colnames(betas) %in% as.vector(SURF1_data$basename)]

# Order betas in same order as Lifespan Data
filtered_betas <- filtered_betas[,match(filtered_data$basename, colnames(filtered_betas))]
Control_betas <- Control_betas[,match(Control_data$basename, colnames(Control_betas))]
SURF1_betas <- SURF1_betas[,match(SURF1_data$basename, colnames(SURF1_betas))]

# Order betas in same order as Lifespan Data
filtered_betas <- filtered_betas[,match(filtered_data$basename, colnames(filtered_betas))]
Control_betas <- Control_betas[,match(Control_data$basename, colnames(Control_betas))]
SURF1_betas <- SURF1_betas[,match(SURF1_data$basename, colnames(SURF1_betas))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_betas)==filtered_data$basename)
identical(colnames(filtered_betas),as.character(filtered_data$basename))


# Load P values from Mixed Effects Model
file_location <- "/Users/gabrielsturm/NYSPI G-drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/SURF1_Analysis/LMER"
setwd(file_location)
dir()
pvalues <- read.csv("LMER_DNAm_SURF1_75_days.csv")

# Load bed file
file_location2 <- "/Users/gabrielsturm/NYSPI G-drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/SURF1_Analysis/DMR"
setwd(file_location2)
dir()

# Make my own Bed File for combp
# A data frame from bed format file with colname name "V1","V2", "V3","V4","V5", 
# V1 indicate chromosome (1,2,3,...,X,Y), V2 is chromosome position, V4 is for P value and V5 for name of CpGs

V1 <- pvalues$chr 
unique(V1)
V1 <- substring(V1,4, 5)
unique(V1)
V2 <- as.numeric(pvalues$pos)
V3 <- as.numeric(pvalues$pos)
#pvalues <- pvalues[pvalues$X %in% pvalues$CpG,] 
#identical(pvalues$X,IlluminaAnnot.EPIC$Name)
V4 <- pvalues$P_values 
V5 <- as.character(pvalues$CpG)
bed_data <- data.frame(V1,V2,V3,V4,V5)
head(bed_data)

# Filter for only sites with p < 0.05
bed_data <- bed_data[bed_data$V4 < 0.99,]

unique(bed_data$V1)
any(is.null(bed_data$V4))
nrow(bed_data)

cores <- detectCores() - 4                   ## Number of cores in your computer
cores
  
DMR_results <- combp(bed_data, dist.cutoff=1000, bin.size=310, seed=0.01,
      region_plot=F, mht_plot=F, nCores=cores)
setwd(file_location2)
dir()
combp_results <- read.csv("resu_combp.csv")
nrow(combp_results)

# Obtain number of CpGs, Plots. Genes and Methylation differences for each DMR
DMR.output <- combp_results
ndmrs <- nrow(DMR.output)

library(foreach)
library(doParallel)
library(parallel)

numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl)

inputs <- 1:ndmrs
dmr_diff <- function(i) {
  gene_list <- vector(mode="character",length=5)
  names(gene_list) <- c("gene1", "gene2", "Methylation_diff", "Methylation_direction", "nCpGs")
  
  
  dmr_chrom <- paste0("chr", DMR.output$chr[i])
  start_pos <- DMR.output$start[i]
  end_pos <- DMR.output$end[i]
  dmr_annotation <- pvalues[pvalues$chr == dmr_chrom,]
  dmr_annotation <- dmr_annotation[dmr_annotation$pos >= start_pos,]
  dmr_annotation <- dmr_annotation[dmr_annotation$pos <= end_pos,]
  # Number of CpGs on the DMR
  nDMRsites <- nrow(dmr_annotation)
  nDMRsites
  dmr_gene <- as.character(unique(dmr_annotation$Gene_name))
  dmr_gene <- dmr_gene[dmr_gene !=""]
  
  # Genes in the DMR
  dmr_gene
  if(!is.na(dmr_gene[1])) {
    if(length(dmr_gene) > 1) {
      gene_list$gene1 <- dmr_gene[1]
      gene_list$gene2 <- dmr_gene[2]
    }
    else {
      gene_list$gene1 <- dmr_gene  
    }
  }  
  
  dmr_betas <- filtered_betas[rownames(filtered_betas) %in% dmr_annotation$CpG,]
  if(is.matrix(dmr_betas) == FALSE) {
    gene_list$nCpGs <- 1
    dmr_Control_betas <- dmr_betas[names(dmr_betas) %in% as.vector(Control_data$basename)]
    dmr_SURF1_betas <- dmr_betas[names(dmr_betas) %in% as.vector(SURF1_data$basename)]
    gene_list$Methylation_diff <- median(dmr_SURF1_betas) - median(dmr_Control_betas)
    if(gene_list$Methylation_diff > 0) {
      gene_list$Methylation_direction <- "Hypermethylated"
    }
    else if(gene_list$Methylation_diff < 0) {
      gene_list$Methylation_direction <- "Hypomethylated"
    }
  }
  else {
    dmr_Control_betas <- dmr_betas[,colnames(dmr_betas) %in% as.vector(Control_data$basename)]
    dmr_SURF1_betas <- dmr_betas[,colnames(dmr_betas) %in% as.vector(SURF1_data$basename)]
    gene_list$Methylation_diff <- median(rowMedians(dmr_SURF1_betas)) - median(rowMedians(dmr_Control_betas))
    if(gene_list$Methylation_diff > 0) {
      gene_list$Methylation_direction <- "Hypermethylated"
    }
    else if(gene_list$Methylation_diff < 0) {
      gene_list$Methylation_direction <- "Hypomethylated"
    }
    gene_list$nCpGs <- nrow(dmr_betas)
  }
  return(gene_list)
}
results <- foreach(i=inputs, .packages='matrixStats') %dopar% {
  dmr_diff(i)
}

on.exit(stopCluster(cl))

gene_list <- do.call(rbind, results)
head(gene_list)
dim(gene_list)


combp_results <- data.frame(combp_results, gene_list)
combp_results$Methylation_diff <- as.numeric(combp_results$Methylation_diff)
combp_results$nCpGs <- as.numeric(combp_results$nCpGs)
combp_results$gene1 <- as.character(combp_results$gene1)
combp_results$gene2 <- as.character(combp_results$gene2)
combp_results$Methylation_direction <- as.character(combp_results$Methylation_direction)
#combp_results <- combp_results[,1:10]
head(combp_results)
dim(combp_results)
nrow(combp_results)
print(hist(combp_results$nCpGs))






# Get negative log p value for every DMR
signed_negative_log_pvalue <- -log(combp_results$p)
for(i in 1:nrow(combp_results)) {
  if(combp_results$Methylation_direction[i] == 'Hypomethylated') {
    signed_negative_log_pvalue[i] <- signed_negative_log_pvalue[i] * -1
  }
}
head(signed_negative_log_pvalue)
combp_results <- data.frame(combp_results, signed_negative_log_pvalue)
head(combp_results)

cutoff <- 10
clusters <- vector(mode ="numeric", length(nrow(combp_results)))
for(i in 1:nrow(combp_results)) {
  value <- as.numeric(combp_results$signed_negative_log_pvalue[i])
  if(value > cutoff) {
    clusters[i] <- 2
    #print(value)
  }
  
  else if(value < -cutoff) {
    clusters[i] <- 0
    #print(value)
  }
  else {
    clusters[i] <- 1
  }
}
length(clusters)
length(clusters[clusters == 2])
length(clusters[clusters == 1])
length(clusters[clusters == 0])
combp_results <- cbind(combp_results, clusters)
head(combp_results)
dim(combp_results)




# Filter for just DMRs with > 2 CpGs
combp_results_filtered <- combp_results[!is.na(combp_results$nCpGs),]
combp_results_filtered <- combp_results_filtered[as.numeric(combp_results_filtered$nCpGs) > 2,]
nrow(combp_results_filtered) # 15,231

setwd(file_location2)
write.csv(combp_results_filtered, "DMRs_Combp_SURF1_75_days_3+dmrs.csv")
write.csv(combp_results, "DMRs_Combp_SURF1_75_days_all_dmrs.csv")

combp_results_filtered <- read.csv("DMRs_Combp_SURF1_75_days_3+dmrs.csv")

# iPAGE input file
iPAGE_input <- data.frame(combp_results_filtered$gene1, combp_results_filtered$signed_negative_log_pvalue)
iPAGE_input <- iPAGE_input[iPAGE_input[,1] !="",]
colnames(iPAGE_input) <- c("","values")
head(iPAGE_input)
nrow(iPAGE_input)
length(unique(iPAGE_input[,1]))
#add back missing genes with no sig DMR
setwd(genewizDir)
Genelist <- as.vector(read.csv("unique_gene_names.csv", row.names = 1))
length(Genelist[,1])
missing_genes <- Genelist[!Genelist[,1] %in% iPAGE_input[,1],]
length(missing_genes)
missing_genes <- data.frame(missing_genes, rep("0", length(missing_genes)))
colnames(missing_genes) <- c("","values")
iPAGE_input <- rbind(iPAGE_input, missing_genes)
head(iPAGE_input)
nrow(iPAGE_input)
setwd("~/Documents/Balaji/DNAm/inputs/")
write.table(iPAGE_input,"DMR_LMER_DNAm_SURF1_75_days.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

iPAGE_cluster_input <- data.frame(combp_results_filtered$gene1, combp_results_filtered$clusters)
iPAGE_cluster_input <- iPAGE_cluster_input[iPAGE_cluster_input[,1] !="",]
colnames(iPAGE_cluster_input) <- c("","values")
head(iPAGE_cluster_input)
#add back missing genes with no sig DMR
setwd(genewizDir)
Genelist <- as.vector(read.csv("unique_gene_names.csv", row.names = 1))
length(Genelist[,1])
missing_genes <- Genelist[!Genelist[,1] %in% iPAGE_cluster_input[,1],]
length(missing_genes)
missing_genes <- data.frame(missing_genes, rep("1", length(missing_genes)))
colnames(missing_genes) <- c("","values")
iPAGE_cluster_input <- rbind(iPAGE_cluster_input, missing_genes)
head(iPAGE_cluster_input)
nrow(iPAGE_cluster_input)
setwd("~/Documents/Balaji/DNAm/inputs/")
write.table(iPAGE_cluster_input,"DMR_LMER_DNAm_SURF1_75_days_clusters.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

file_location2 <- "/Users/gabrielsturm/NYSPI G-drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/DNA Methylation/Part 2/SURF1_Analysis/DMR"
setwd(file_location2)
combp_results_filtered <- read.csv("DMRs_Combp_SURF1_75_days_3+dmrs.csv", row.names = 1)
combp_results <- read.csv("DMRs_Combp_SURF1_75_days_all_dmrs.csv", row.names = 1)
head(combp_results_filtered)

DMR.output <- combp_results_filtered
ndmrs <- nrow(DMR.output)
ndmrs
pdf(file=paste0(file_location2,"/DMRs_SURF1_75_days_for_publication.pdf"), width=2.2, height=1.5) 
for(i in 1:10) {
  # Get all CpGs in Top DMR
  #i <- 87
  dmr_chrom <- paste0("chr", DMR.output$chr[i])
  start_pos <- DMR.output$start[i]
  end_pos <- DMR.output$end[i]
  dmr_annotation <- IlluminaAnnot.EPIC[IlluminaAnnot.EPIC$chr == dmr_chrom,]
  dmr_annotation <- dmr_annotation[dmr_annotation$pos >= start_pos,]
  dmr_annotation <- dmr_annotation[dmr_annotation$pos <= end_pos,]
  # Number of CpGs on the DMR
  nDMRsites <- nrow(dmr_annotation)
  nDMRsites
  dmr_gene <- as.character(unique(dmr_annotation$corrected_gene_name))
  
  # Genes in the DMR
  dmr_gene
  dmr_betas <- filtered_betas[rownames(filtered_betas) %in% dmr_annotation$X,]
  if(is.matrix(dmr_betas) == FALSE) {
    next;
  }
  rownames(dmr_betas)
  dmr_annotation$X
  dmr_data_1 <- data.frame(dmr_betas, dmr_annotation)
  library(gdata)
  unbeta <- unmatrix(dmr_betas, byrow=TRUE)
  names(unbeta)
  # Repeat filtered data for every CpG Site
  repeated_filtered_data <- do.call("rbind", replicate(nDMRsites, filtered_data[,1:19], simplify = FALSE))
  nSamples <- nrow(filtered_data)
  #rep.row(dmr_annotation, nSamples)
  repeated_dmr_annotation <-  dmr_annotation[rep(seq_len(nrow(dmr_annotation)), each = nSamples), ]
  dim(repeated_dmr_annotation)
  dmr_data <- data.frame(repeated_filtered_data, names(unbeta), unbeta, repeated_dmr_annotation)
  title <- paste0("Gene: ", unique(as.character(repeated_dmr_annotation$corrected_gene_name)),
                  ", ", unique(as.character(repeated_dmr_annotation$chr)), ", ", (DMR.output$end[i] - DMR.output$start[i]), "bp")
  #Site_info <-  paste(dmr_annotation$X, dmr_annotation$corrected_gene_group, dmr_annotation$Relation_to_Island, sep=", ")
  Site_info <-  paste(substring(dmr_annotation$corrected_gene_group,0,2), substring(dmr_annotation$Relation_to_Island,0,2), sep="_")
  basesize=7
  text_size <- basesize
  if(nDMRsites  > 12) {
    text_size <- basesize/nDMRsites *basesize
    if(text_size<3) {
      text_size <- 3
    }
  }
  DMR_plot <- ggplot(data=dmr_data, aes(y=unbeta,x=X, group = Clinical_Condition, color = Clinical_Condition), fill = Clinical_Condition) +
    geom_point(alpha =0.4) + #, aes(shape =Cell_Line)) +
    #geom_line(position = 'jitter') +
    #geom_smooth(method="loess") +
    stat_summary(fun=mean, geom="line", aes(color = Clinical_Condition), size = 1, alpha =1) +
    scale_fill_manual(values=c("darkgray", "blue4")) +
    scale_color_manual(values=c("darkgray", "blue4")) +
    scale_y_continuous(name="Beta", limits=c(0,1)) +
    #geom_area(mapping = aes(y = ifelse(unbeta>0.1 & unbeta< 0.4 , unbeta, 0)), fill = "red") +
    scale_x_discrete(name="", labels= Site_info) +
    ggtitle(title) +
    #annotate("text", label = ), x=1,  y = 0.5, size = 8) +
    theme_classic(base_size = basesize) +
    theme(axis.text.x = element_text(size = text_size,vjust = 1.01,hjust=1,angle = 55), #  , 
          legend.position="none",
          axis.line = element_line(colour = 'black', size = 0.3),
          #axis.text=element_text(size=8),
          legend.title = element_blank())
  print(DMR_plot)
}
dev.off() 

#group <- c(rep(10,1),rep(12,5),rep(14,20),rep(16,26),rep(18,15),rep(20,32),rep(22,24),rep(24,27),rep(26,4),rep(28,2),rep(30,1))
#sd(group) / mean(group) * 100

# Volcano Plot
# SURF1 Volcano Plot
library(ggrepel)
volcano_plot <- ggplot(combp_results, aes(x = Methylation_diff, y = -log(p))) + 
  geom_point(aes(color=Methylation_direction), alpha = 0.3) +
  scale_color_manual(values = c("red", "blue","gray")) +
  geom_hline(yintercept=-log(0.01), linetype="dashed", color ="black", size = 0.5) +
  geom_vline(xintercept=0, size=0.5) +
  scale_x_continuous(name = "Methylation Difference") +
  scale_y_continuous(name = "-log(p value)") +
  geom_text_repel(data=subset(combp_results_filtered, abs(Methylation_diff) > .7 | -log(p) > 80),
                  box.padding = unit(0.1, "lines"),alpha = 0.5, segment.alpha = 0.5, #min.segment.length = Inf,
                  aes(x = Methylation_diff, y = -log(p), label=gene1), hjust=0, vjust=0) +
  theme_classic(base_size = 20) +
  theme(legend.title = element_blank(), 
        legend.position="bottom", 
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.text=element_text(size=6))
volcano_plot
setwd(file_location2)
pdf("Volcano_plot_DMR_DNAm_SURF1_75_days.pdf")
volcano_plot
# Close the pdf file
dev.off()


# Proportion Upregulated vs Downregulated
downR <- nrow(combp_results_filtered[combp_results_filtered$Methylation_direction == "Hypomethylated",])
downR # 8,141
upR <- nrow(combp_results_filtered[combp_results_filtered$Methylation_direction == "Hypermethylated",])
upR # 7,090
unR <- nrow(combp_results_filtered[combp_results_filtered$Methylation_direction == "Unchanged",])
unR
downR_prop <- round(downR / (downR+upR+unR) * 100,1)
downR_prop # 54%
upR_prop <- round(upR / (downR+upR+unR) * 100,1)
upR_prop # 47%
unR_prop <- round(unR / (downR+upR+unR) * 100,1)
unR_prop # 0%


### Pi Chart ###
library(ggplot2)
library(scales)
library(tidyverse)
# Barplot
group <- c("Hypomethylated", "Hypermethylated")
value <- c(downR, upR)
percentC <- c(downR_prop, upR_prop)
df <- data.frame(group, value, percentC)
df

# Compute the position of labels
df <- df %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(df$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pie <- ggplot(df, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  
  geom_text(aes(y = ypos, label = paste0(group, "\n", percentC,"%")), color = "white", size=8) +
  scale_fill_manual(values=c("red","blue"))
pie

setwd(file_location2)
pdf("Pichart_DMR_DNAm_SURF1_75_days.pdf")
pie
dev.off()