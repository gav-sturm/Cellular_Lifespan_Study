library('Rtsne')
library('viridis')
library('ggplot2')
library('grid')
library('gridExtra')
library(devtools)
library('Vennerable') # install_github("js229/Vennerable")

# source('~/Dropbox/Columbia/scripts/enrichGO_Symbol.R')
# annot = read.csv('~/Dropbox/RNAseq_Lifespan_Sample_info_ed.txt', sep='\t', header=T)
# annot2 = annot; rownames(annot2) = paste(annot$Sample.ID,'_lenscaleTPM',sep='')
# annot2$sample_name = rownames(annot2)

# oligo_meta = read.csv('/data2/mito_data/postprocess/kallisto/data_iPAGE_fromGS/iPAGE_data_for_balaji/Oligo_sample_metadata.csv', sep=',', header=T)
# oligo_samps = data.frame(Sample.ID = as.character(unique(oligo_meta$RNAseq_ID)))
# surf_meta = read.csv('/data2/mito_data/postprocess/kallisto/data_iPAGE_fromGS/iPAGE_data_for_balaji/SURF1_sample_metadata.csv', sep=',', header=T)
# surf_samps = data.frame(Sample.ID = as.character(unique(surf_meta$RNAseq_ID)))


# oligo_lmer = read.csv('/data2/mito_data/postprocess/kallisto/data_iPAGE_fromGS/iPAGE_data_for_balaji/Oligo_LMER_results_transcripts.csv', sep=',', header=T)
# 
# surf_lmer = read.csv('/data2/mito_data/postprocess/kallisto/data_iPAGE_fromGS/iPAGE_data_for_balaji/SURF1_LMER_results_genes.csv', sep=',', header=T)

# oligo_meta = unique(merge(annot2, oligo_samps)) ; surf_meta = unique(merge(annot2, surf_samps))
# rownames(oligo_meta) = as.character(oligo_meta$sample_name)
# rownames(surf_meta) = as.character(surf_meta$sample_name)

genewizDir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/RNA/RNAseq/Genewiz/"
setwd(genewizDir)
dir()

files <- list.files('kallisto_output_newv3/')
files <- files[-352]
files
f_ = paste('kallisto_output_newv3/', files,'/abundance.tsv',sep='')
names(f_) = paste(files,'_lenscaleTPM',sep='')

x = read.csv('kallisto_output_newv3/1/abundance.tsv',sep='\t',header=T)
g_ = sapply(lapply(strsplit(as.character(x$target_id), '\\;'), 'rev'), '[', 2)
tx2gene = data.frame(target_id = as.character(x$target_id), gene = g_)

# https://www.rdocumentation.org/packages/tximport/versions/1.0.3/topics/tximport
txi = tximport(f_, type='kallisto', tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM")
txi2 = tximport(f_, type='kallisto', tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM", txOut=TRUE)

gene_count_df <- txi$counts
gene_count_df[1:5,1:5]
gene_abundance_df <- txi$abundance
gene_abundance_df[1:5,1:5]
gene_length_df <- txi$length
gene_length_df[1:5,1:5]
setwd(genewizDir)
write.csv(gene_count_df,"Balaji_DEseq/txi_gene_count_lengthScaledTPM.csv")
write.csv(gene_abundance_df,"Balaji_DEseq/txi_gene_abundance.csv")
write.csv(gene_length_df,"Balaji_DEseq/txi_gene_length.csv")

transcript_count_df <- txi2$counts
transcript_count_df[1:5,1:5]
transcript_abundance_df <- txi2$abundance
transcript_abundance_df[1:5,1:5]
transcript_length_df <- txi2$length
transcript_length_df[1:5,1:5]
write.csv(transcript_count_df,"Balaji_DEseq/txi_transcript_lengthScaledTPM.csv")
write.csv(transcript_abundance_df,"Balaji_DEseq/txi_transcript_abundance.csv")
write.csv(transcript_length_df,"Balaji_DEseq/txi_transcript_length.csv")

setwd(genewizDir)
gene_count_df <- read.csv("Balaji_DEseq/txi_gene_count_lengthScaledTPM.csv")
colnames(gene_count_df) <- substring(colnames(gene_count_df),2)

setwd(genewizDir)
annot2 <- read.csv('Balaji_DEseq/sample_annotation', sep='\t', header=T)
annot_df = merge(annot2, data.frame(sample_name = colnames(gene_count_df))) ; rownames(annot_df) = as.character(annot_df$sample_name)

# save annotation
setwd(genewizDir)
write.csv(annot_df,"Balaji_DEseq/Balaji_RNAseq_annotation.csv")

# LS annotation
ls_metadata <- LS_Data[!is.na(LS_Data$RNAseq_ID),]
ls_metadata <-ls_metadata[order(ls_metadata$RNAseq_ID),]
nrow(ls_metadata)
setwd(genewizDir)
write.csv(ls_metadata,"Balaji_DEseq/GEO_RNAseq_annotation2.csv")

setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/data_sharing/")
LS_metadata <- read.csv("Cellular_lifespan_study_metadata.csv",row.names = 1,header=T)
LS_metadata <-LS_metadata[!is.na(LS_metadata$RNAseq_sampleID),]
LS_metadata <-LS_metadata[order(LS_metadata$RNAseq_sampleID),]
nrow(LS_metadata)
setwd(genewizDir)
write.csv(LS_metadata,"Balaji_DEseq/GEO_RNAseq_annotation2.csv")



# samples removed by QC
removed_samples <- seq(1:360)[!seq(1:360) %in% annot_df$Selected_Samples]
removed_samples
#samples excluded after sequencing
excluded_samples <- annot_df$Selected_Samples[!annot_df$Selected_Samples %in% LS_metadata$RNAseq_sampleID]  
excluded_samples

removed_samples <- data.frame(removed_samples, rep("did not pass pre-seq QC",length(removed_samples)))
colnames(removed_samples) <- c("sample","reason")
excluded_samples <- data.frame(excluded_samples, rep("excluded post-seq",length(excluded_samples)))
colnames(excluded_samples) <- c("sample","reason")
missing_samples <- rbind(removed_samples, excluded_samples)
missing_samples

# missing_samples metadata
all_sample_metadata <-read.csv("RNAseq_Lifespan_Sample_info.csv")
dim(all_sample_metadata)
missing_samples_metadata <- all_sample_metadata[all_sample_metadata$Sample.Name %in% missing_samples$sample,]
nrow(missing_samples_metadata)
missing_samples_metadata <- missing_samples_metadata[match(missing_samples$sample,missing_samples_metadata$Sample.Name),]
identical(missing_samples$sample, missing_samples_metadata$Sample.Name)
missing_samples <- data.frame(missing_samples, missing_samples_metadata)

write.csv(missing_samples,"Balaji_DESeq/excluded_samples_RNAseq.csv")

# metadata for GEO sumbission
setwd("/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Gabriel Sturm/Projects/Project 2- Cell Lifespan Aging/shinyapp")
dir()
GEO_metadata <- read.csv("downloadable_data/Cellular_Lifespan_study_RNAseq.csv")
GEO_metadata <- GEO_metadata[!is.na(GEO_metadata$RNAseq_sample),]
GEO_metadata <- GEO_metadata[order(GEO_metadata$RNAseq_sample),]
rownames(GEO_metadata) <- GEO_metadata$RNAseq_sampleID
nrow(GEO_metadata)

setwd(genewizDir)
write.csv(GEO_metadata,"GEO_RNAseq_metadata.csv")

# final processed matrix for GEO submission
setwd(genewizDir)
vst_matrix <- as.matrix(read.table("vst_gene_level_GS",header=T))
head(vst_matrix)
max(vst_matrix)
min(vst_matrix)
dim(vst_matrix)

hist(as.matrix(vst_matrix), breaks=100)
# Rename Samples
old_sample_names <- colnames(vst_matrix)
corrected_sample_name <- substring(colnames(vst_matrix),2,str_locate(colnames(vst_matrix), "_")[,1]-1)
corrected_sample_name <- paste0("Sample_",corrected_sample_name)
corrected_sample_name
colnames(vst_matrix) <- corrected_sample_name
head(vst_matrix)
max(vst_matrix)
min(vst_matrix)
dim(vst_matrix)

# Remove excluded samples
final_vst_matrix <- vst_matrix[,!colnames(vst_matrix) %in% paste0("Sample_",missing_samples$sample)]
dim(final_vst_matrix)

# put samples in order
geo_dir <- "/Users/gabrielsturm/NYSPI G-Drive/MitoLab - General/ Members Folders/Sharing files/For Gabriel/SURF1 mito disease paper/GEO_submission"
setwd(geo_dir)
GEO_metadata <- read.csv("RNAseq/GEO_RNAseq_metadata.csv")
final_vst_matrix <- final_vst_matrix[,match(GEO_metadata$Sample.name,colnames(final_vst_matrix))]
write.csv(final_vst_matrix,"processed_cell_lifespan_RNAseq_data.csv")

# Raw counts matrix (kallisto output)
setwd(genewizDir)
raw_counts_df <- as.matrix(read.table("Balaji_DEseq/est_counts_kallisto",header=T,row.names=1))
# rownames(raw_counts_df) <- raw_counts_df[,1]
# raw_counts_df <- raw_counts_df[,-1]
# raw_counts_df <- as.matrix(sapply(raw_counts_df, as.numeric))
raw_counts_df[1:5,1:5]
dim(raw_counts_df)
max(raw_counts_df)
min(raw_counts_df)

# Rename Samples
old_sample_names <- colnames(raw_counts_df)
corrected_sample_name <- substring(colnames(raw_counts_df),2,str_locate(colnames(raw_counts_df), "_")[,1]-1)
corrected_sample_name <- paste0("Sample_",corrected_sample_name)
corrected_sample_name
colnames(raw_counts_df) <- corrected_sample_name
head(raw_counts_df)
max(raw_counts_df)
min(raw_counts_df)
dim(raw_counts_df)

# Remove excluded samples
final_raw_counts_matrix <- raw_counts_df[,!colnames(raw_counts_df) %in% paste0("Sample_",missing_samples$sample)]
dim(final_raw_counts_matrix)

# order samples
setwd(geo_dir)
final_raw_counts_matrix <- final_raw_counts_matrix[,match(GEO_metadata$Sample.name,colnames(final_raw_counts_matrix))]
head(final_raw_counts_matrix)
max(final_raw_counts_matrix)
min(final_raw_counts_matrix)
dim(final_raw_counts_matrix)
write.csv(final_raw_counts_matrix,"raw_counts_cell_lifespan_RNAseq_data.csv")

# get MD5 checksums
setwd(genewizDir)
md5s <- read.table("Balaji_DEseq/30-359114334.md5")
# sort them
file_names <- substring(md5s$V2,str_locate(md5s$V2, "fastq/")[,1]+6)
sample_names <- paste0("Sample_",substring(file_names,1,str_locate(file_names, "_")[,1]-1))
md5s <- data.frame(sample_names, md5s, file_names)
dim(md5s)
# remove extra samples
# Remove excluded samples
md5s <- md5s[!md5s$sample_names %in% paste0("Sample_",missing_samples$sample),]
dim(md5s)
# remove files not in metadata
md5s <- md5s[md5s$sample_names %in% GEO_metadata$Sample.name,]
dim(md5s)
# reorder samples
file_order <- c(GEO_metadata$raw.file,str_replace(GEO_metadata$raw.file, "R1", "R2"))
length(file_order)
md5s <- md5s[match(file_order,md5s$file_names),]
dim(md5s)
write.csv(md5s, "Balaji_DEseq/fastq_md5s.csv")



# source('~/Dropbox/Columbia/scripts/enrichGO_Symbol.R')
# annot = read.csv('~/Dropbox/RNAseq_Lifespan_Sample_info_ed.txt', sep='\t', header=T)
# annot2 = annot; rownames(annot2) = paste(annot$Sample.ID,'_lenscaleTPM',sep='')
# annot2$sample_name = rownames(annot2)
#load lifespan study data



dds_new = DESeqDataSetFromTximport(txi, colData=annot_df, design= ~ Clinical.Condition + Sex + Treatment + Donor.Age + Passage)
#dds_new_DE = DESeq(dds_new)
#vst_dds_new_DE = vst(dds_new)
vst_df_full = as.data.frame(assay(vst(dds_new)))

#f_ = paste('/data2/mito_data/kallisto_output_newv3/', list.files('/data2/mito_data/kallisto_output_newv3/'),'/abundance.tsv',sep='')
#names(f_) = paste(list.files('/data2/mito_data/kallisto_output_newv3/'),'_lenscaleTPM',sep='')
f_surf = f_[intersect(as.character(surf_meta$sample_name), names(f_))]
f_oligo = f_[intersect(as.character(oligo_meta$sample_name), names(f_))]
f_surf_oligo = f_[intersect(union(as.character(surf_meta$sample_name), as.character(oligo_meta$sample_name)), names(f_))]
s_o_meta = merge(data.frame(sample_name = names(f_surf_oligo)), annot2)
rownames(s_o_meta) = as.character(s_o_meta$sample_name)

txi_surf = tximport(f_surf, type='kallisto', tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM")
txi_oligo = tximport(f_oligo, type='kallisto', tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM")
txi_s_o = tximport(f_surf_oligo, type='kallisto', tx2gene=tx2gene, countsFromAbundance = "lengthScaledTPM")

dds_surf = DESeqDataSetFromTximport(txi_surf, colData=surf_meta, design= ~ Clinical.Condition + Time.Grown..days.)
dds_oligo = DESeqDataSetFromTximport(txi_oligo, colData=oligo_meta, design= ~ Treatment + Cell.Line + Time.Grown..days.)
dds_s_o = DESeqDataSetFromTximport(txi_s_o, colData=s_o_meta[colnames(txi_s_o$counts),], design= ~ Clinical.Condition + Treatment + Time.Grown..days.)


surf_DE = DESeq(dds_surf) ; oligo_DE = DESeq(dds_oligo)
vst_s_o = vst(dds_s_o)
vst_df = as.data.frame(assay(vst_s_o))

mds_ = as.data.frame(cmdscale(as.dist(1-cor(vst_df, method='spearman')))) ; mds_$sample_name = rownames(mds_)
mds_ = merge(mds_, s_o_meta)

pdf('~/Dropbox/mito_data/tsne_mito.pdf', width=15, height=10)
for(perp_ in c(5:13)) {
  t_ = Rtsne(unique(t(vst_df)),perplexity=perp_)
  ts_ = as.data.frame(t_$Y) ; ts_$sample_name = colnames(vst_df) ; ts_ = merge(ts_, s_o_meta)
  p1 = ggplot(ts_, aes(V1, V2, label=Cell.Line, shape=paste(Treatment,Clinical.Condition))) + geom_point(aes(color=Time.Grown..days.),size=4) +scale_colour_gradient(low='blue', high='red') + ggtitle(paste('tSNE', perp_))
  p2 = ggplot(ts_, aes(V1, V2, label=Cell.Line, shape=paste(Treatment,Clinical.Condition))) + geom_point(aes(color=Time.Grown..days.),size=4) +scale_colour_gradient(low='blue', high='red') + geom_text() + ggtitle(paste('tSNE',perp_))
  grid.arrange(p1, p2, nullGrob(), nullGrob())
}
ggplot(mds_, aes(V1, V2, label=Cell.Line, shape=paste(Treatment,Clinical.Condition))) + geom_point(aes(color=Time.Grown..days.),size=4) +scale_colour_gradient(low='blue', high='red') + ggtitle('MDS')
ggplot(mds_, aes(V1, V2, label=Cell.Line, shape=paste(Treatment,Clinical.Condition))) + geom_point(aes(color=Time.Grown..days.),size=4) +scale_colour_gradient(low='blue', high='red') + geom_text() + ggtitle('MDS')

dev.off()

res_surf = as.data.frame(results(surf_DE)) ; res_oligo = as.data.frame(results(oligo_DE))
res_surf$Gene = rownames(res_surf) ; res_oligo$Gene = rownames(res_oligo)

surf_rep = res_surf[which(res_surf$log2FoldChange < 0 & res_surf$padj < 0.05),]
surf_act = res_surf[which(res_surf$log2FoldChange > 0 & res_surf$padj < 0.05),]
oligo_rep = res_oligo[which(res_oligo$log2FoldChange < 0 & res_oligo$padj < 0.05),]
oligo_act = res_oligo[which(res_oligo$log2FoldChange > 0 & res_oligo$padj < 0.05),]

l_surf = list(SURF1_rep = as.character(surf_rep$Gene), 
              SURF1_act = as.character(surf_act$Gene),
              SURF_nde = setdiff(as.character(res_surf$Gene), union(as.character(surf_rep$Gene), as.character(surf_act$Gene))))
l_oligo = list(Oligo_rep = as.character(oligo_rep$Gene), 
               Oligo_act = as.character(oligo_act$Gene),
               Oligo_nde = setdiff(as.character(res_oligo$Gene), union(as.character(oligo_rep$Gene), as.character(oligo_act$Gene))))
l_s_o1 = list(SURF1_rep = as.character(surf_rep$Gene), 
              SURF1_act = as.character(surf_act$Gene),
              SURF_nde = setdiff(as.character(res_surf$Gene), union(as.character(surf_rep$Gene), as.character(surf_act$Gene))),
              Oligo_rep = as.character(oligo_rep$Gene), 
              Oligo_act = as.character(oligo_act$Gene),
              Oligo_nde = setdiff(as.character(res_oligo$Gene), union(as.character(oligo_rep$Gene), as.character(oligo_act$Gene))))

l_s_o2 = list(SURF1_rep = as.character(surf_rep$Gene), 
              SURF1_act = as.character(surf_act$Gene),
              Oligo_rep = as.character(oligo_rep$Gene), 
              Oligo_act = as.character(oligo_act$Gene))

s_act = as.character(surf_act$Gene) ; s_not_act = setdiff(as.character(res_surf$Gene), s_act)
o_act = as.character(oligo_act$Gene) ; o_not_act = setdiff(as.character(res_oligo$Gene), o_act)
s_rep = as.character(surf_rep$Gene) ; s_not_rep = setdiff(as.character(res_surf$Gene), s_rep)
o_rep = as.character(oligo_rep$Gene) ; o_not_rep = setdiff(as.character(res_oligo$Gene), o_rep)
m_act = rbind(c(length(intersect(s_act, o_act)), length(intersect(s_act, o_not_act))),
              c(length(intersect(s_not_act, o_act)), length(intersect(s_not_act, o_not_act))))
m_rep = rbind(c(length(intersect(s_rep, o_rep)), length(intersect(s_rep, o_not_rep))),
              c(length(intersect(s_not_rep, o_rep)), length(intersect(s_not_rep, o_not_rep))))
fish_act = fisher.test(m_act) ; fish_rep = fisher.test(m_rep)
if (fish_act$p.value < 10^-125) { p_a = '<1e-100'}
if (fish_rep$p.value < 10^-125) { p_r = '<1e-100'}

pdf('~/Dropbox/mito_data/pie_Venn_SURF1_Oligo_DE.pdf', width=10, height=10)
plot(Venn(l_s_o1), doWeights=F) ; plot(Venn(l_s_o2), doWeights=F, type='squares')
pie(as.numeric(lapply(l_surf, length)), labels=paste(names(l_surf), '(n=', as.numeric(lapply(l_surf, length)),')',sep=''), main='SURF1 DESeq2')
pie(as.numeric(lapply(l_oligo, length)), labels=paste(names(l_oligo), '(n=', as.numeric(lapply(l_oligo, length)),')',sep=''), main='Oligo DESeq2')
dev.off()



ipage_surf = data.frame(Gene = as.character(res_surf$Gene))
ipage_surf$categ = 1
ipage_surf$categ[which(res_surf$log2FoldChange < 0 & res_surf$padj < 0.05)] = 2
ipage_surf$categ[which(res_surf$log2FoldChange > 0 & res_surf$padj < 0.05)] = 0

ipage_oligo = data.frame(Gene = as.character(res_oligo$Gene))
ipage_oligo$categ = 1
ipage_oligo$categ[which(res_oligo$log2FoldChange < 0 & res_oligo$padj < 0.05)] = 2
ipage_oligo$categ[which(res_oligo$log2FoldChange > 0 & res_oligo$padj < 0.05)] = 0

write.table(ipage_surf, '/data2/mito_data/postprocess/kallisto/DESeq_SURF_2021-05-07', sep='\t', quote=F, row.names=F)
write.table(ipage_oligo, '/data2/mito_data/postprocess/kallisto/DESeq_Oligo_2021-05-07', sep='\t', quote=F, row.names=F)


surf_go_up = enrichGO_Symbol(data.frame(Approved.Symbol = as.character(surf_act$Gene))); surf_go_do = enrichGO_Symbol(data.frame(Approved.Symbol = as.character(surf_rep$Gene)))
surf_go_u = surf_go_up[[2]] ; surf_go_d = surf_go_do[[2]]
surf_go_u = surf_go_u[sort(surf_go_u$classic, index.return=T)$ix,] ; surf_go_d = surf_go_d[sort(surf_go_d$classic, index.return=T)$ix,]

oligo_go_up = enrichGO_Symbol(data.frame(Approved.Symbol = as.character(oligo_act$Gene))); oligo_go_do = enrichGO_Symbol(data.frame(Approved.Symbol = as.character(oligo_rep$Gene)))
oligo_go_u = oligo_go_up[[2]] ; oligo_go_d = oligo_go_do[[2]]
oligo_go_u = oligo_go_u[sort(oligo_go_u$classic, index.return=T)$ix,] ; oligo_go_d = oligo_go_d[sort(oligo_go_d$classic, index.return=T)$ix,]


#export PAGEDIR=/home/balaji/tavazoie_lab/software/iPAGE/
#perl /home/balaji/tavazoie_lab/software/iPAGE/page.pl --expfile=DESeq_SURF_2021-05-07 --species=human_gene_symbols --exptype=discrete --independence=0 --catmincount=20
#perl /home/balaji/tavazoie_lab/software/iPAGE/page.pl --expfile=DESeq_Oligo_2021-05-07 --species=human_gene_symbols --exptype=discrete --independence=0 --catmincount=20
