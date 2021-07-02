
library(limma) 
library(ggplot2)
library(mgcv)
library(lme4)
library(tidyverse)
library(ggrepel)

genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
dir()

vst_matrix <- as.matrix(read.table("vst_gene_level_GS",header=T))
head(vst_matrix)
max(vst_matrix)
min(vst_matrix)
dim(vst_matrix)
hist(as.matrix(vst_matrix), breaks=100)
# Rename Samples
old_sample_names <- colnames(vst_matrix)
corrected_sample_name <- substring(colnames(vst_matrix),2,str_locate(colnames(vst_matrix), "_")[,1]-1)
colnames(vst_matrix) <- corrected_sample_name
Gene_data <- 2^vst_matrix
min(Gene_data)
max(Gene_data)


Gene_names <- rownames(Gene_data)
head(Gene_names)
length(Gene_names)


# Add 1 to all tpm values to avoid infinite fold change
Transcript_data <- log2(Gene_data)
max(rowMax(Transcript_data))
min(rowMin(Transcript_data))

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
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 75,]

# Split data by control vs SURF1 samples
Control_data <- filtered_data[filtered_data$Clinical_Condition == "Normal",]
SURF1_data <- filtered_data[filtered_data$Clinical_Condition == "SURF1_Mutation",]

# Timepoints
Control_timepoints <- Control_data$Days_Grown
SURF1_timepoints <- SURF1_data$Days_Grown

# Retrieve TPM for control and SURF1 samples
# Genes
# filtered_genes <- Gene_data[,colnames(Gene_data) %in% as.vector(filtered_data$RNAseq_ID)]
# Control_genes <- Gene_data[,colnames(Gene_data) %in% as.vector(Control_data$RNAseq_ID)]
# SURF1_genes <- Gene_data[,colnames(Gene_data) %in% as.vector(SURF1_data$RNAseq_ID)]
# Transcripts
filtered_transcripts <- Transcript_data[,colnames(Transcript_data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_transcripts <- Transcript_data[,colnames(Transcript_data) %in% as.vector(Control_data$RNAseq_ID)]
SURF1_transcripts <- Transcript_data[,colnames(Transcript_data) %in% as.vector(SURF1_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
# Genes
# filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
# Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
# SURF1_genes <- SURF1_genes[,match(SURF1_data$RNAseq_ID, colnames(SURF1_genes))]
# Transcripts
filtered_transcripts <- filtered_transcripts[,match(filtered_data$RNAseq_ID, colnames(filtered_transcripts))]
Control_transcripts <- Control_transcripts[,match(Control_data$RNAseq_ID, colnames(Control_transcripts))]
SURF1_transcripts <- SURF1_transcripts[,match(SURF1_data$RNAseq_ID, colnames(SURF1_transcripts))]



## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_transcripts)==filtered_data$RNAseq_ID)
identical(colnames(filtered_transcripts),as.character(filtered_data$RNAseq_ID))
# all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
# identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))


library(foreach)
library(doParallel)
library(parallel)

numCores <- detectCores() - 2
cl <- makeCluster(numCores)
registerDoParallel(cl)

inputs <- 1:nrow(filtered_transcripts)
LMER_test <- function(i) {
  gene_tpm <- as.numeric(filtered_transcripts[i,])
  #Test Gene
  #gene_tpm <- as.numeric(filtered_genes[59,])
  P_values <- 1
  # if(sum(gene_tpm) < 15) {
  #   P_values <- 1
  # }
  #else {
    Gene <- data.frame(filtered_data, gene_tpm)
    
    MEmodel = lmer(gene_tpm ~ Clinical_Condition*Days_Grown + (1|Days_Grown),
                   data = Gene,  REML = F)
    summary <- summary(MEmodel)
    # Gene %>%
    #   # save predicted values
    #   mutate(pred_dist = fitted(MEmodel)) %>%
    #   # graph
    #   ggplot(aes(x=Days_Grown, y=pred_dist, group=Clinical_Condition, color=Clinical_Condition)) + theme_classic() +
    #   geom_line(size=1)
    
    
    null_model <- lmer(gene_tpm ~ Days_Grown + (1|Days_Grown),
                       data = Gene,  REML = F)
    
    anova <- anova(null_model, MEmodel)
    P_values <- anova$`Pr(>Chisq)`[2]
  #}
  
  
  message(i)
  
  #Organize matrix
  results <- data.frame(rownames(filtered_transcripts)[i], P_values)
  return(results)
}
results <- foreach(i=inputs, .packages='lme4') %dopar% {
  LMER_test(i)
}

stopCluster(cl)

#mresults <- as.data.frame(matrix(unlist(results), nrow=nrow(filtered_betas), ncol=2))
lmer_results <- do.call(rbind, results)
#mresults <- as.data.frame(t(as.data.frame(results)))
head(lmer_results)
results <- lmer_results

# Correct for multipule Comparisons
bonfi_P_values <- as.matrix(p.adjust(results[,2], method = "bonferroni", n = length(results[,2])))
fdr_P_values <- as.matrix(p.adjust(results[,2], method = "fdr", n = length(results[,2])))

#Organize matrix
results <- data.frame(Gene_names, results, bonfi_P_values, fdr_P_values)
colnames(results)[1] <- "Gene"
#colnames(results)[3] <- "FullGeneName"
#head(results)
dim(results)

#setwd(genewizDir)
#write.csv(results, "LMER/LMER_allRNA_SURF1_75_days_transcripts_1tpmcutoff_logtransformed.csv")
#results <- read.csv("LMER/LMER_allRNA_SURF1_75_days_transcripts_1tpmcutoff_logtransformed.csv")

# Function for caluating the methylation differnce between two lists of betas
expression_difference <- function(group1_betas, group2_betas) {
  # Split Sites by hyper and hypo methylation SURF1 vs control
  Methyl_diff <- matrix(ncol=3, nrow=nrow(group1_betas))
  for(i in 1:nrow(group1_betas)) {
    diff <- median(group2_betas[i,]) - median(group1_betas[i,])
    Methyl_diff[i,2] <- diff
    if(diff > 0) {
      Methyl_diff[i,1] <- "Upregulated"
      Methyl_diff[i,3] <- -log(results$P_values[i]) * 1
    }
    else if (diff < 0) {
      Methyl_diff[i,1] <- "Downregulated"
      Methyl_diff[i,3] <- -log(results$P_values[i]) * -1
    }
    else {
      Methyl_diff[i,1] <- "Unchanged"
      Methyl_diff[i,3] <- -log(results$P_values[i]) * 1
    }
    #Methyl_diff[i,3] <- median(group2_betas[i,]) / median(group1_betas[i,])
    
  }
  rownames(Methyl_diff) <- rownames(group1_betas)
  colnames(Methyl_diff) <- c("expression_direction", "expression_fold_change", "signed_negative_log_pvalue")
  return(Methyl_diff)
}

# Run methyl differnce comparing control vs SURF1
expression_diff <- expression_difference(Control_transcripts, SURF1_transcripts) # 1min
head(expression_diff)
expression_diff <-  expression_diff[order(match(rownames(expression_diff),results$FullGeneName)),]  
identical(rownames(expression_diff), results$Gene)
results <- data.frame(results, expression_diff)
#results$expression_median_difference <- as.numeric(results$expression_median_difference)
results$expression_fold_change <- as.numeric(results$expression_fold_change)
results$signed_negative_log_pvalue <- as.numeric(results$signed_negative_log_pvalue)
# make all unchanged genes to a pvalue of 1
results[results$expression_direction == "Unchanged",3] <- 1
results[results$expression_direction == "Unchanged",4] <- 1
results[results$expression_direction == "Unchanged",5] <- 1
results[results$expression_direction == "Unchanged",8] <- 1
table(results$expression_direction)

head(results)
max(results$signed_negative_log_pvalue)
min(results$signed_negative_log_pvalue)
hist(results$signed_negative_log_pvalue, breaks = 50)

# Make 3 manual clusters of upregulated, downregulated, and unchanged genes
cutoff <- 2
clusters <- vector(mode ="numeric", length(nrow(results)))
for(i in 1:nrow(results)) {
  value <- as.numeric(results$signed_negative_log_pvalue[i])
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
results <- cbind(results, clusters)
head(results)
dim(results)

### Plot the top sites of results
#order sites by category
results <- results[order(results$P_values),]
head(results)

# choose number of sites to plot
plot_top_sites <- function(n_sites, results, linear) {
  top_sites <- results$Gene[1:n_sites]
  for(i in 1:length(top_sites)) {
    Gene <- top_sites[i]
    #Gene <- "ELOVL2"
    Site_betas <- c(Control_transcripts[rownames(Control_transcripts) == Gene,],SURF1_transcripts[rownames(SURF1_transcripts) == Gene,]) 
    Site_Data <- cbind(filtered_data,Site_betas)
    
    p <- ggplot(Site_Data, aes(x=Days_Grown, y = Site_betas)) +
      geom_point(size = 2 , stroke = 0.5, aes(shape = Cell_Line, fill = Cell_Line)) +
      #geom_line(aes(y = Control_GAM)) +
      scale_color_manual(values = c("dimgray", "blue")) +
      scale_fill_manual(values =  c("dimgray","dimgray","dimgray","blue","blue","blue")) +
      scale_shape_manual(values=c(21,22,23,21,22,23)) +
      scale_y_continuous(name = Gene) +
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
      control_gam <- gam(Control_transcripts[rownames(Transcript_data) == Gene,]~s(Control_timepoints,bs="cr",k=3))
      SURF1_gam <- gam(SURF1_transcripts[rownames(Transcript_data) == Gene,]~s(SURF1_timepoints,bs="cr",k=3)) 
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
plot_top_sites(10, results, T)

setwd(genewizDir)
write.csv(results, "LMER/LMER_allRNA_SURF1_75_days_vst.csv")

setwd(genewizDir)
results <- read.csv("LMER/LMER_allRNA_SURF1_75_days_vst.csv")


iPAGE_input <- data.frame(results$Gene, results$signed_negative_log_pvalue)
colnames(iPAGE_input) <- c("","values")
setwd("~/Documents/Balaji/inputs/")
write.table(iPAGE_input,"LMER_allRNA_SURF1_75_days_vst.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

iPAGE_cluster_input <- data.frame(results$Gene, results$clusters)
colnames(iPAGE_cluster_input) <- c("","values")
write.table(iPAGE_cluster_input,"LMER_allRNA_SURF1_75_days_genes_vst_clusters.txt", append = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)


# Attach results to annotation
cutoff <- 0.05
sig_results <- results[results$fdr_P_values < cutoff,]
nrow(sig_results) # 7,142

setwd(genewizDir)
write.csv(sig_results, "LMER/LMER_allRNA_SURF1_75_days_vst_sig_genes.csv")


# Volcano Plot
# SURF1 Volcano Plot
library(ggrepel)
volcano_plot <- ggplot(results, aes(x = expression_fold_change, y = -log(P_values))) + 
  geom_point(aes(color=expression_direction), alpha = 0.3) +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_hline(yintercept=-log(0.05), linetype="dashed", color ="black", size = 0.5) +
  geom_vline(xintercept=0, size=0.5) +
  scale_x_continuous(name = "Log 2 TPM Fold Change") +
  scale_y_continuous(name = "-log(p value)") +
  geom_text_repel(data=subset(sig_results, abs(expression_fold_change) > 3.5 | -log(P_values) > 25), max.overlaps =100,
                  box.padding = unit(0.1, "lines"),alpha = 0.5, segment.alpha = 0.5, #min.segment.length = Inf,
                  aes(x = expression_fold_change, y = -log(P_values), label=Gene), hjust=0, vjust=0) +
  theme_classic(base_size = 20) +
  theme(legend.title = element_blank(), 
        legend.position="bottom", 
        axis.line = element_line(colour = 'black', size = 0.5),
        axis.text=element_text(size=6))
volcano_plot
setwd(genewizDir)
pdf("LMER/Volcano_plot_allRNA_SURF1_75_days_vst.pdf")
volcano_plot
# Close the pdf file
dev.off() 


# Proportion Upregulated vs Downregulated
downR <- nrow(sig_results[sig_results$expression_direction == "Downregulated",])
downR
upR <- nrow(sig_results[sig_results$expression_direction == "Upregulated",])
upR
unR <- nrow(sig_results[sig_results$expression_direction == "Unchanged",])
unR
downR_prop <- round(downR / (downR+upR+unR) * 100,1)
downR_prop # 81%
upR_prop <- round(upR / (downR+upR+unR) * 100,1)
upR_prop # 19%
unR_prop <- round(unR / (downR+upR+unR) * 100,1)
unR_prop # 0%


### Pi Chart ###
library(ggplot2)
library(scales)
# Barplot
group <- c("Downregulated", "Upregulated")
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
  scale_fill_manual(values=c("blue","red"))
pie

setwd(genewizDir)
pdf("LMER/Pichart_allRNA_SURF1_75_days_vst.pdf")
pie
dev.off()