library(ggplot2)
library(tidyverse)
library(qlcMatrix)
library(matrixStats)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
dir()

# Gene_data <- as.matrix(read.csv("allRNA_Gene_RNAseq_data_no_cutoff.csv", row.names = 1))
# # Rename Samples
# old_sample_names <- colnames(Gene_data)
# corrected_sample_name <- substring(colnames(data),5,7)
# colnames(data) <- corrected_sample_name
# # Add 1 to all tpm values to avoid infinite fold change
# Gene_data <- Gene_data + 1
# max(rowMax(Gene_data))
# min(rowMin(Gene_data))

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

# Get Gene Lists
setwd(genewizDir)
gene_lists <- t(read.csv("RNAseq_Gene_Lists.csv", row.names =1))


#load lifespan study data
baseDir3 <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp"
setwd(baseDir3)
dir()
LS_Data <- read.csv("Lifespan_Study_Data.csv")

data <- Gene_data

### Heatmap Function ###
graphHeatmap <- function(heatmap_data, genes_tpm_data, gene_group, group_names, group_lines, group_size, group_colors, ColSeps) {
  nGroups <- length(group_names)
  colors <- brewer.pal(n = 5, name = "Spectral")
  colors <- factor(c("#000099", "#3399FF" ,"#FAFAFA", "#FF8000", "#FF0000"), levels = c("#000099", "#3399FF" ,"#FAFAFA", "#FF8000", "#FF0000"))
  
  # Map Colors to Max Min of heatmap data with same magnitude in both directions
  mapping <- max(heatmap_data)
  if(abs(min(heatmap_data)) > max(heatmap_data)) {
    mapping <- abs(min(heatmap_data))
  }
  if(mapping > 4) {
    mapping <- 4
  }
  col_fun = colorRamp2(c(-mapping,-mapping/2,0,mapping/2,mapping), colors)
  
  Groups <- vector()
  for(i in 1:nGroups) {
    Groups <-  c(Groups, rep(group_names[i], group_size[i]))
  }
  Groups <- factor(Groups, levels = group_names)
  #message(Groups)
  colnames(heatmap_data) <- Groups
  
  Groups_col <- vector()
  for(i in 1:nGroups) {
    Groups_col <-  c(Groups_col, rep(group_colors[i], group_size[i]))
  }
  Groups_col <- factor(Groups_col, levels = group_colors)
  #message(Groups_col)
  
  Groups2 <- factor(filtered_data$Cell_Line_Group,levels = unique(filtered_data$Cell_Line_Group))
  nrep <- 3
  if(unique(filtered_data$Study_Part) == 4) {
    nrep <- 2
  }
  Groups3 <- vector()
  for(i in 1:nGroups) {
    Groups3 <-  c(Groups3, rep(group_names[i], nrep))
  }
  Groups3 <- factor(Groups3, levels = group_names)
  
  Groups3_col <- vector()
  for(i in 1:nGroups) {
    Groups3_col <-  c(Groups3_col, rep(group_colors[i], nrep))
  }
  Groups3_col <- as.vector(factor(Groups3_col, levels = group_colors))
  
  Cell_Lines <- group_lines
  row_txt_size <- 10 - (nrow(heatmap_data) *.1)
  if(row_txt_size < 4) {
    row_txt_size <- 4
  }
  label_size <- 8 - (nGroups * 1)
  # lowex_genes <- genes_tpm_data[rowMedians(genes_tpm_data) < 1,]
  # lowex_at <- which(rowMedians(genes_tpm_data) < 1)
  # update(lowex_at)
  # lowex_labels <- rep('low expression', nrow(lowex_genes))
  # update(lowex_labels)
  heatmap <- Heatmap(heatmap_data, col=col_fun, cluster_rows = TRUE, cluster_columns = FALSE,
                     column_title = gene_group,
                     row_names_gp = gpar(fontsize = row_txt_size),
                     column_title_gp = gpar(fontsize = 16, fontface = "bold"), 
                     cluster_column_slices = FALSE, row_names_side = "left", 
                     show_column_names = FALSE, #column_title = NULL,
                     top_annotation = HeatmapAnnotation(
                       Cell_Lines = anno_block(labels= Cell_Lines,labels_gp = gpar(fontsize = label_size), gp = gpar(fill = "white")),
                       Days_Grown = anno_barplot(filtered_data$Days_Grown),
                       show_annotation_name = T,
                       #annotation_name = c(rep("Days \n Grown", nGroups)),
                       annotation_name_side = "left", 
                       annotation_name_gp = gpar(fontsize = 6)
                       #annotation_name_offset = unit(0, "mm")
                     ),
                     right_annotation = rowAnnotation(
                       TPM = anno_boxplot(genes_tpm_data,width = unit(1, "cm"), outline = FALSE),
                       show_annotation_name = T
                       # if(lowex_labels > 0) {
                       #   low_expression = anno_mark(at = lowex_at, labels = lowex_labels)
                       # }
                     ),
                     bottom_annotation = HeatmapAnnotation(
                       Treatment = anno_block(labels = Groups3, labels_gp = gpar(fontsize = label_size), gp = gpar(fill = Groups3_col))),
                     name = "log2 mcTPM", column_split = Groups2, column_gap = unit(ColSeps, "mm"))
  return(heatmap)
}

#### GeneSet Heatmap Function ###
GeneSet_Heatmaps <- function() {
  npathways_missing <- 0
  pathways_missing <- vector()
  for(i in 1:ncol(gene_lists)) {
    genes <- gene_lists[,i]
    genes <- genes[genes != ""]
    if(length(genes) > 1 & !grepl("http", genes[1]) &  !grepl("\\s+", genes[1])) {
      genes_data <- Log2_Median_Centered_Data[rownames(Log2_Median_Centered_Data) %in% genes,]
      genes_tpm_data <- filtered_genes[rownames(filtered_genes) %in% genes,]
      nMissing <- length(genes[!genes %in% rownames(genes_data)])
      print(paste0(colnames(gene_lists)[i]," # of Missing Genes: ", nMissing))
      if(dim(genes_data)[1] > 0) {
        heatmap <- graphHeatmap(genes_data, genes_tpm_data, colnames(gene_lists)[i], group_names, group_lines, group_size, group_colors, ColSeps)
        print(heatmap)
      }
      else {
        pathways_missing <- c(pathways_missing, colnames(gene_lists[,i]))
        npathways_missing <- npathways_missing +1
      }
    }
  }
  print(paste0("# of Missing Pathways: ", npathways_missing))
  print(paste0("Missing Pathways: ", paste(pathways_missing,collapse=", ")))
}


### Breakup by Experiment ###

##### SURF1 & Oligo #######

# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14","hFB6","hFB7","hFB8")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatments %in% c("Control", "Oligomycin"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs SURF1 vs Oligo samples
Control_data <- filtered_data[filtered_data$Clinical_Condition == "Normal",]
Control_data <- Control_data[Control_data$Treatments == "Control",]
nrow(Control_data)
SURF1_data <- filtered_data[filtered_data$Clinical_Condition == "SURF1_Mutation",]
SURF1_data <- SURF1_data[SURF1_data$Treatments == "Control",]
nrow(SURF1_data)
Oligo_data <- filtered_data[filtered_data$Treatments == "Oligomycin",]
Oligo_data <- Oligo_data[Oligo_data$Clinical_Condition == "Normal",]
nrow(Oligo_data)

Control_timepoints <- Control_data$Days_Grown
SURF1_timepoints <- SURF1_data$Days_Grown
Oligo_timepoints <- Oligo_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
SURF1_genes <- data[,colnames(data) %in% as.vector(SURF1_data$RNAseq_ID)]
Oligo_genes <- data[,colnames(data) %in% as.vector(Oligo_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
SURF1_genes <- SURF1_genes[,match(SURF1_data$RNAseq_ID, colnames(SURF1_genes))]
Oligo_genes <- Oligo_genes[,match(Oligo_data$RNAseq_ID, colnames(Oligo_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

### Median Center all Genes to the youngest samples of the controls
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0

# Get the median of the youngest control sample
young_average <- Control_genes[,1:3]
for(j in 1:length(unique(Control_data$Cell_Line))) {
  cell_line <- unique(Control_data$Cell_Line)[j]
  line_samples <- Control_data[Control_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- Control_genes[,colnames(Control_genes) %in% line_samples$RNAseq_ID]
  
  young_average[,j] <- line_genes[,1]
  start <- end + 1
}
young_average <- rowMedians(young_average)

# center each cell line to its youngest control and SURF1 to the median of all the controls
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  if(cell_line == "hFB6"| cell_line == "hFB7" | cell_line == "hFB8") {
    young_centered <- line_genes / young_average
  }
  else {
    young_centered <- line_genes / line_genes[,1]
  }
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}
head(Median_Centered_Data)



# Center all samples to median of control samples
# Median_Centered_Data <- matrix(ncol=ncol(filtered_genes), nrow=nrow(filtered_genes))
# for(i in 1:nrow(filtered_genes)) {
#   Median_Centered_Data[i,] <- filtered_genes[i,] / median(Control_genes[i,])
# }
# rownames(Median_Centered_Data) <- rownames(filtered_genes)
# colnames(Median_Centered_Data) <- colnames(filtered_genes)
# head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)

group_names <- c("Control", "SURF1", "Oligo")
group_lines <- c("Donor1","Donor2","Donor3","Patient1","Patient2","Patient3","Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(SURF1_data), nrow(Oligo_data))
group_colors <- c("gray", "royalblue1", "darkorchid1")
ColSeps <- c(rep(0.8,2),3,rep(0.8,2),3,rep(0.8,2))


## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps/Heatmaps_SURF1_Oligo_RNAseq_3_vst.pdf")
GeneSet_Heatmaps()
dev.off()

###############

##### Aging Controls ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatments == "Control",]
filtered_data <-  filtered_data[filtered_data$Study_Part %in% c(1,2),]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs SURF1 vs Oligo samples
filtered_data <- filtered_data[filtered_data$Clinical_Condition == "Normal",]
#Control_data <- Control_data[Control_data$Treatments == "Control",]
nrow(filtered_data)

#Control_timepoints <- Control_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
#Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
#Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))
#write.csv(cbind(filtered_data, filtered_genes[rownames(filtered_genes) == 'YME1L1']), 'YME1L1.csv')
# Add 1 to all tpm values to avoid infinite fold change
filtered_genes <- filtered_genes + 1
max(rowMax(filtered_genes))
min(rowMin(filtered_genes))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)

max(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)

group_names <- c("Control")
group_lines <- c("Donor1","Donor2","Donor3")
group_size <- c(nrow(filtered_data))
group_colors <- c("azure3")
ColSeps <- c(rep(0.8,2))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Aging_Controls_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################


##### Contact Inhibition ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "Contact_Inhibition_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part %in% c(2,3),]
filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatments == "Control",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatments == "Contact_Inhibition",]
nrow(CI_data)

Control_timepoints <- Control_data$Days_Grown
CI_timepoints <- CI_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))


# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control", "Contact \n Inhibition")
group_lines <- c("Donor2","Donor3","Donor4", "Donor2","Donor3","Donor4")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "dodgerblue1")
ColSeps <- c(rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Contact_Inhibition_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Hypoxia ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "Control_3"),]
filtered_data <-  filtered_data[filtered_data$Study_Part ==3 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "Control_3",]
nrow(CI_data)

Control_timepoints <- Control_data$Days_Grown
CI_timepoints <- CI_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))


# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control", "Hypoxia")
group_lines <- c("Donor2","Donor3","Donor4", "Donor2","Donor3","Donor4")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "deeppink3")
ColSeps <- c(rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Hypoxia_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Contact Inhibition + Hypoxia ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB11","hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21","Control_3","Contact_Inhibition_21", "Contact_Inhibition_3"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 3 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 120,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)
Hypoxia_data <- filtered_data[filtered_data$Treatment == "Control_3",]
nrow(Hypoxia_data)
CI_data <- filtered_data[filtered_data$Treatment == "Contact_Inhibition_21",]
nrow(CI_data)
CI_Hypoxia_data <- filtered_data[filtered_data$Treatment == "Contact_Inhibition_3",]
nrow(CI_Hypoxia_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
Hypoxia_genes <- data[,colnames(data) %in% as.vector(Hypoxia_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]
CI_Hypoxia_genes <- data[,colnames(data) %in% as.vector(CI_Hypoxia_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
Hypoxia_genes <- Hypoxia_genes[,match(Hypoxia_data$RNAseq_ID, colnames(Hypoxia_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]
CI_Hypoxia_genes <- CI_Hypoxia_genes[,match(CI_Hypoxia_data$RNAseq_ID, colnames(CI_Hypoxia_genes))]


## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))


# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control", "Hypoxia","Contact \n Inhibition","Contact \n Inhibition \n Hypoxia")
group_lines <- c("Donor2","Donor3","Donor4", "Donor2","Donor3","Donor4",
                 "Donor2","Donor3","Donor4", "Donor2","Donor3","Donor4")
group_size <- c(nrow(Control_data), nrow(Hypoxia_data), nrow(CI_data), nrow(CI_Hypoxia_data))
group_colors <- c("azure3", "deeppink3","blue","green")
ColSeps <- c(rep(0.8,2),5,rep(0.8,2),5,rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Contact_Inhibition_Hypoxia_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### DEX ######
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "DEX_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "DEX_21",]
nrow(CI_data)

Control_timepoints <- Control_data$Days_Grown
CI_timepoints <- CI_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# save data for Natalia
colnames(filtered_genes) <- filtered_data$Unique_Variable_Name
head(filtered_genes)
# setwd(genewizDir)
# write.csv(filtered_genes, "DEX_RNAseq_Gene_Data_nocutoff.csv")

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "DEX")
group_lines <- c("Donor1","Donor2","Donor3", "Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "red")
ColSeps <- c(rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_DEX_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### SURF1 + DEX ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14","hFB6","hFB7","hFB8")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21","DEX_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]
nrow(filtered_data)

# Split data by control vs surf1 vs dex vs surf1 + dex
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
Control_data <- Control_data[Control_data$Clinical_Condition == "Normal",]
nrow(Control_data)
SURF1_data <- filtered_data[filtered_data$Treatment == "Control_21",]
SURF1_data <- SURF1_data[SURF1_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_data)
DEX_data <- filtered_data[filtered_data$Treatment == "DEX_21",]
DEX_data <- DEX_data[DEX_data$Clinical_Condition == "Normal",]
nrow(DEX_data)
SURF1_DEX_data <- filtered_data[filtered_data$Treatment == "DEX_21",]
SURF1_DEX_data <- SURF1_DEX_data[SURF1_DEX_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_DEX_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
SURF1_genes <- data[,colnames(data) %in% as.vector(SURF1_data$RNAseq_ID)]
DEX_genes <- data[,colnames(data) %in% as.vector(DEX_data$RNAseq_ID)]
SURF1_DEX_genes <- data[,colnames(data) %in% as.vector(SURF1_DEX_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
ncol(filtered_genes)
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
SURF1_genes <- SURF1_genes[,match(SURF1_data$RNAseq_ID, colnames(SURF1_genes))]
DEX_genes <- DEX_genes[,match(DEX_data$RNAseq_ID, colnames(DEX_genes))]
SURF1_DEX_genes <- SURF1_DEX_genes[,match(SURF1_DEX_data$RNAseq_ID, colnames(SURF1_DEX_genes))]


## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))


# Get the median of the youngest control sample
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0

young_average <- Control_genes[,1:3]
for(j in 1:length(unique(Control_data$Cell_Line))) {
  cell_line <- unique(Control_data$Cell_Line)[j]
  line_samples <- Control_data[Control_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- Control_genes[,colnames(Control_genes) %in% line_samples$RNAseq_ID]
  
  young_average[,j] <- line_genes[,1]
  start <- end + 1
}
young_average <- rowMedians(young_average)

# center each cell line to its youngest control and SURF1 to the median of all the controls
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  #print(line_samples$Cell_Line_Group)
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  if(cell_line == "hFB6"| cell_line == "hFB7" | cell_line == "hFB8") {
    young_centered <- line_genes / young_average
  }
  else {
    young_centered <- line_genes / line_genes[,1]
  }
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}
head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control","SURF1","DEX","SURF1 \n DEX")
group_lines <- c("Donor1","Donor2","Donor3","Patient1","Patient2","Patient3",
                 "Donor1","Donor2","Donor3","Patient1","Patient2","Patient3")
group_size <- c(nrow(Control_data), nrow(SURF1_data), nrow(DEX_data), nrow(SURF1_DEX_data))
group_colors <- c("azure3", "royalblue1", "red", "pink")
ColSeps <- c(rep(0.8,2),5,rep(0.8,2),5,rep(0.8,2),5,rep(0.8,3))
ncol(Log2_Median_Centered_Data)
## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_SURF1_DEX_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Modulators ######
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "Modulators_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "Modulators_21",]
nrow(CI_data)

Control_timepoints <- Control_data$Days_Grown
CI_timepoints <- CI_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "Modulators")
group_lines <- c("Donor1","Donor2","Donor3", "Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "green")
ColSeps <- c(rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Modulators_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Modulators + DEX ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21","Modulators_21", "DEX_21", "Modulators+DEX_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]
nrow(filtered_data)

# Split data by control vs surf1 vs dex vs surf1 + dex
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
#Control_data <- Control_data[Control_data$Clinical_Condition == "Normal",]
nrow(Control_data)
SURF1_data <- filtered_data[filtered_data$Treatment == "Modulators_21",]
#SURF1_data <- SURF1_data[SURF1_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_data)
DEX_data <- filtered_data[filtered_data$Treatment == "DEX_21",]
DEX_data <- DEX_data[DEX_data$Clinical_Condition == "Normal",]
nrow(DEX_data)
SURF1_DEX_data <- filtered_data[filtered_data$Treatment == "Modulators+DEX_21",]
#SURF1_DEX_data <- SURF1_DEX_data[SURF1_DEX_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_DEX_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
SURF1_genes <- data[,colnames(data) %in% as.vector(SURF1_data$RNAseq_ID)]
DEX_genes <- data[,colnames(data) %in% as.vector(DEX_data$RNAseq_ID)]
SURF1_DEX_genes <- data[,colnames(data) %in% as.vector(SURF1_DEX_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
ncol(filtered_genes)
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
SURF1_genes <- SURF1_genes[,match(SURF1_data$RNAseq_ID, colnames(SURF1_genes))]
DEX_genes <- DEX_genes[,match(DEX_data$RNAseq_ID, colnames(DEX_genes))]
SURF1_DEX_genes <- SURF1_DEX_genes[,match(SURF1_DEX_data$RNAseq_ID, colnames(SURF1_DEX_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Get the median of the youngest control sample
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0

young_average <- Control_genes[,1:3]
for(j in 1:length(unique(Control_data$Cell_Line))) {
  cell_line <- unique(Control_data$Cell_Line)[j]
  line_samples <- Control_data[Control_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- Control_genes[,colnames(Control_genes) %in% line_samples$RNAseq_ID]
  
  young_average[,j] <- line_genes[,1]
  start <- end + 1
}
young_average <- rowMedians(young_average)

# center each cell line to its youngest control and SURF1 to the median of all the controls
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  #print(line_samples$Cell_Line_Group)
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  if(cell_line == "hFB6"| cell_line == "hFB7" | cell_line == "hFB8") {
    young_centered <- line_genes / young_average
  }
  else {
    young_centered <- line_genes / line_genes[,1]
  }
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}
head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control","Modulators","DEX","Modualators \n DEX")
group_lines <- c("Donor1","Donor2","Donor3","Donor1","Donor2","Donor3",
                 "Donor1","Donor2","Donor3","Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(SURF1_data), nrow(DEX_data), nrow(SURF1_DEX_data))
group_colors <- c("azure3", "green", "red", "darkgreen")
ColSeps <- c(rep(0.8,2),5,rep(0.8,2),5,rep(0.8,2),5,rep(0.8,3))
ncol(Log2_Median_Centered_Data)
## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Modulators_DEX_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Oligomycin ######
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "Oligomycin_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "Oligomycin_21",]
nrow(CI_data)

Control_timepoints <- Control_data$Days_Grown
CI_timepoints <- CI_data$Days_Grown

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "Oligomycin")
group_lines <- c("Donor1","Donor2","Donor3", "Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "darkorchid1")
ColSeps <- c(rep(0.8,2),5,rep(0.8,3))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Oligomycin_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Oligo + DEX ######
# Filter Lifespan Data to hFB12,13,14 and hFB6,7,8 Control
selected_cell_lines <- c("hFB12","hFB13","hFB14")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21","Oligomycin_21", "DEX_21", "Oligomycin+DEX_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 2 ,]
#filtered_data <-  filtered_data[filtered_data$Replicate_Line != 3,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]
nrow(filtered_data)

# Split data by control vs surf1 vs dex vs surf1 + dex
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
#Control_data <- Control_data[Control_data$Clinical_Condition == "Normal",]
nrow(Control_data)
SURF1_data <- filtered_data[filtered_data$Treatment == "Oligomycin_21",]
#SURF1_data <- SURF1_data[SURF1_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_data)
DEX_data <- filtered_data[filtered_data$Treatment == "DEX_21",]
DEX_data <- DEX_data[DEX_data$Clinical_Condition == "Normal",]
nrow(DEX_data)
SURF1_DEX_data <- filtered_data[filtered_data$Treatment == "Oligomycin+DEX_21",]
#SURF1_DEX_data <- SURF1_DEX_data[SURF1_DEX_data$Clinical_Condition == "SURF1_Mutation",]
nrow(SURF1_DEX_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
SURF1_genes <- data[,colnames(data) %in% as.vector(SURF1_data$RNAseq_ID)]
DEX_genes <- data[,colnames(data) %in% as.vector(DEX_data$RNAseq_ID)]
SURF1_DEX_genes <- data[,colnames(data) %in% as.vector(SURF1_DEX_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
ncol(filtered_genes)
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
SURF1_genes <- SURF1_genes[,match(SURF1_data$RNAseq_ID, colnames(SURF1_genes))]
DEX_genes <- DEX_genes[,match(DEX_data$RNAseq_ID, colnames(DEX_genes))]
SURF1_DEX_genes <- SURF1_DEX_genes[,match(SURF1_DEX_data$RNAseq_ID, colnames(SURF1_DEX_genes))]

## Check data columns on M-values align with Phenotype as in SentrixID_Position on Chip
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Get the median of the youngest control sample
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0

young_average <- Control_genes[,1:3]
for(j in 1:length(unique(Control_data$Cell_Line))) {
  cell_line <- unique(Control_data$Cell_Line)[j]
  line_samples <- Control_data[Control_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- Control_genes[,colnames(Control_genes) %in% line_samples$RNAseq_ID]
  
  young_average[,j] <- line_genes[,1]
  start <- end + 1
}
young_average <- rowMedians(young_average)

# center each cell line to its youngest control and SURF1 to the median of all the controls
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  #print(line_samples$Cell_Line_Group)
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  if(cell_line == "hFB6"| cell_line == "hFB7" | cell_line == "hFB8") {
    young_centered <- line_genes / young_average
  }
  else {
    young_centered <- line_genes / line_genes[,1]
  }
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}
head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)

group_names <- c("Control","Oligomycin","DEX","Oligomycin \n DEX")
group_lines <- c("Donor1","Donor2","Donor3","Donor1","Donor2","Donor3",
                 "Donor1","Donor2","Donor3","Donor1","Donor2","Donor3")
group_size <- c(nrow(Control_data), nrow(SURF1_data), nrow(DEX_data), nrow(SURF1_DEX_data))
group_colors <- c("azure3", "darkorchid1", "red", "darkorchid3")
ColSeps <- c(rep(0.8,2),5,rep(0.8,2),5,rep(0.8,2),5,rep(0.8,3))
ncol(Log2_Median_Centered_Data)
## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Oligomycin_DEX_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################


##### 2-Deoxyglucose ######
selected_cell_lines <- c("hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "2-Deoxyglucose_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 4 ,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "2-Deoxyglucose_21",]
nrow(CI_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns match
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "2-Deoxyglucose")
group_lines <- c("Donor1","Donor2", "Donor1","Donor2")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "yellow")
ColSeps <- c(rep(0.8,1),5,rep(0.8,2))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_2-Deoxyglucose_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Beta-hydroxybutarate ######
selected_cell_lines <- c("hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "betahydroxybutyrate_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 4 ,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "betahydroxybutyrate_21",]
nrow(CI_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns match
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "Betahydroxybutyrate")
group_lines <- c("Donor1","Donor2", "Donor1","Donor2")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "purple")
ColSeps <- c(rep(0.8,1),5,rep(0.8,2))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Betahydroxybutyrate_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################

##### Galactose_21 ######
selected_cell_lines <- c("hFB12","hFB13")
filtered_data <- LS_Data[LS_Data$Cell_Line %in% selected_cell_lines,]
filtered_data <-  filtered_data[filtered_data$Treatment %in% c("Control_21", "Galactose_21"),]
filtered_data <-  filtered_data[filtered_data$Study_Part == 4 ,]
filtered_data <-  filtered_data[!is.na(filtered_data$RNAseq_ID),]
filtered_data$RNAseq_ID
# Filter for first 50 days of samples of each cell line
#filtered_data <- filtered_data[filtered_data$Days_Grown < 150,]

# Split data by control vs CI vs Oligo samples
Control_data <- filtered_data[filtered_data$Treatment == "Control_21",]
nrow(Control_data)

CI_data <- filtered_data[filtered_data$Treatment == "Galactose_21",]
nrow(CI_data)

# Retrieve TPM for control and Oligo samples
filtered_genes <- data[,colnames(data) %in% as.vector(filtered_data$RNAseq_ID)]
Control_genes <- data[,colnames(data) %in% as.vector(Control_data$RNAseq_ID)]
CI_genes <- data[,colnames(data) %in% as.vector(CI_data$RNAseq_ID)]

# Order betas in same order as Lifespan Data
filtered_genes <- filtered_genes[,match(filtered_data$RNAseq_ID, colnames(filtered_genes))]
Control_genes <- Control_genes[,match(Control_data$RNAseq_ID, colnames(Control_genes))]
CI_genes <- CI_genes[,match(CI_data$RNAseq_ID, colnames(CI_genes))]

## Check data columns match
all(colnames(filtered_genes)==filtered_data$RNAseq_ID)
identical(colnames(filtered_genes),as.character(filtered_data$RNAseq_ID))

# Matricized - Center all genes to the youngest sample of each cell line
Median_Centered_Data <- filtered_genes
start <- 1
end <- 0
for(j in 1:length(unique(filtered_data$Cell_Line))) {
  cell_line <- unique(filtered_data$Cell_Line)[j]
  line_samples <- filtered_data[filtered_data$Cell_Line == cell_line,]
  end <- start + nrow(line_samples) - 1
  line_genes <- filtered_genes[,colnames(filtered_genes) %in% line_samples$RNAseq_ID]
  print(colnames(line_genes))
  young_centered <- line_genes / line_genes[,1]
  
  Median_Centered_Data[,colnames(Median_Centered_Data) %in% line_samples$RNAseq_ID] <- young_centered
  start <- end + 1
}

head(Median_Centered_Data)

Log2_Median_Centered_Data <- log2(Median_Centered_Data)
head(Log2_Median_Centered_Data)
min(Log2_Median_Centered_Data)
max(Log2_Median_Centered_Data)


group_names <- c("Control", "Galactose")
group_lines <- c("Donor1","Donor2", "Donor1","Donor2")
group_size <- c(nrow(Control_data), nrow(CI_data))
group_colors <- c("azure3", "orange")
ColSeps <- c(rep(0.8,1),5,rep(0.8,2))

## Generate Heatmpas for each Gene List
setwd(genewizDir)
pdf("Heatmaps_Galactose_RNAseq_3.pdf")
GeneSet_Heatmaps()
dev.off()
################


