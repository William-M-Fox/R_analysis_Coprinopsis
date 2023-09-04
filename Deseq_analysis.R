#load packages needed
library(viridis)
library(DESeq2)
library(ggdendro)
library(gridExtra)
library(grid)
library(gtable)
library(topGO)
library(GO.db)
library(biomaRt)
library(Rgraphviz)
library(dplyr)
library(tidyr)
library(tibble)
library(gplots)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

#load combined counts data
setwd("/Volumes/student_users/williamfox/R_paper/R_paper/")
counts_stranded <- read.table("counts_stranded.txt", header = TRUE)
counts_reverse <- read.table("counts_stranded_rev.txt", header = TRUE)
counts_stranded_total <- as.data.frame(colSums(counts_stranded[,7:37], na.rm = FALSE))
counts_reverse_total <- as.data.frame(colSums(counts_reverse[,7:37], na.rm = FALSE))
counts_all_total <- cbind(counts_reverse_total, counts_stranded_total)

# DESEQ analysis
fcData <- counts_reverse
rownames(fcData) <- fcData$Geneid
fcData <- fcData[,7:37]
colnames(fcData) <- gsub("^CC_|\\.sam.bam$", "", colnames(fcData))
metaData <- read.csv('CC_IDs.csv')
colnames(metaData) <- c("ID", "stage", "loc", "sample")
metaData$ID <- gsub("^CC_|\\.sam$", "", metaData$ID)
metaData$sample <- gsub("^CC_|\\.sam$", "", metaData$sample)

rownames(metaData) <- metaData$ID

# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData = fcData, 
                              colData = metaData, 
                              design = ~ sample)

# Run DESeq
dds <- DESeq(dds)

# Create a list to store results for each sample
results_list <- list()

for (sample in unique(metaData$sample)) {
  # Set the reference level for this comparison
  dds$sample <- relevel(dds$sample, ref = sample)
  
  # Re-run DESeq
  dds <- DESeq(dds)
  
  # Extract results
  res <- results(dds)
  
  # Remove rows with NA in padj or log2FoldChange
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), ]
  
  # Filter for significant genes
  sig_res_up <- res[res$padj < 0.05 & res$log2FoldChange > 4, ]
  sig_res_down <- res[res$padj < 0.05 & res$log2FoldChange < -4, ]
  
  # Store in the list
  results_list[[sample]] <- list("Up" = sig_res_up, "Down" = sig_res_down)
  combined <- rbind(cbind(sig_res_up, "Direction" = "Up"), cbind(sig_res_down, "Direction" = "Down"))
  
  # Write to a CSV file
  write.csv(combined, file.path("Results/", paste0("deseq_results_", sample, ".csv")), row.names = TRUE)
}

vst_transformer <- vst(dds, blind = FALSE)
vst_data <- assay(vst_transformer)
vst_data <- as.data.frame(vst_data)
colnames(vst_data) <- metaData$ID
vst_data_sub <- vst_data

##PCA plot of vst data
pdf("Results/pca_Plot.pdf", width = 6, height = 6)
pca_plot <- plotPCA(vst_transformer, intgroup = "sample", ntop = 500, returnData = FALSE)
dev.off()

## process vst data for further analysis
long_vst_data_sub <- vst_data_sub %>%
  rownames_to_column("gene_id") %>%
  pivot_longer(cols = -gene_id, names_to = "ID", values_to = "expression")

long_vst_data_sub <- left_join(long_vst_data_sub, metaData, by = "ID")

avg_vst_data_sub <- long_vst_data_sub %>%
  group_by(sample, gene_id) %>%
  summarise(expression = mean(expression, na.rm = TRUE)) %>%
  pivot_wider(names_from = sample, values_from = expression)

avg_vst_data_sub <- as.data.frame(avg_vst_data_sub)
rownames(avg_vst_data_sub) <- avg_vst_data_sub$gene_id
avg_vst_data_sub <- avg_vst_data_sub[, c("VM", "H", "P1", "P2", "YFB_S", "YFB_L", "YFB_C", "FB_S" ,"FB_CL")]

library(purrr)

# Function to process each DESeqResults object
process_DESeqResults <- function(item, name, group) {
  as.data.frame(item) %>% 
    rownames_to_column("Row") %>% 
    mutate(Name = name, Group = group)
}

# Initialize empty dataframes
up_df <- data.frame()
down_df <- data.frame()

# Loop over each element in the up_list and down_list
for (name in names(results_list)) {
  up_df <- bind_rows(up_df, process_DESeqResults(results_list[[name]]$Up, name, "up"))
  down_df <- bind_rows(down_df, process_DESeqResults(results_list[[name]]$Down, name, "down"))
}

#Identfy coding and non coding 4fold DE genes
noncoding_names <- read.csv("Results/noncoding_transcripts.txt", header = FALSE)
noncoding_names$V1 <- gsub("transcript_id ", "", noncoding_names$V1)
noncoding_names <- noncoding_names$V1

coding_names <- read.csv("Results/coding_transcripts.txt", header = FALSE)
coding_names$V1 <- gsub("transcript_id " , "", coding_names$V1)
coding_names <- coding_names$V1

nc_up <- up_df[up_df$Row %in% noncoding_names,]
nc_down <- down_df[down_df$Row %in% noncoding_names,]
c_up <- up_df[up_df$Row %in% coding_names,]
c_down <- down_df[down_df$Row %in% coding_names,]

reg_transcripts <- rbind(nc_up, nc_down, c_up, c_down)

#Vm and primordium specifc
nc_up_VM <- nc_up[nc_up$Name == "VM" |nc_up$Name == "ZM_VM", ]
nc_up_VM_counts <- counts_reverse[counts_reverse$Geneid %in% nc_up_VM$Row,]

nc_up_P1 <- nc_up[nc_up$Name == "P1", ]
nc_up_P1_counts <- counts_reverse[counts_reverse$Geneid %in% nc_up_P1$Row,]

nc_up_P2 <- nc_up[nc_up$Name == "P2", ]
nc_up_P2_counts <- counts_reverse[counts_reverse$Geneid %in% nc_up_P2$Row,]

reg_transcripts_expr <- avg_vst_data_sub[rownames(avg_vst_data_sub) %in% reg_transcripts$Row,]

column_order <- colnames(reg_transcripts_expr)
reg_transcripts_expr <- reg_transcripts_expr[, c(column_order)]

#heatmap of diff expr genes 1
all_diff_genes_heatmap <- heatmap.2(as.matrix(reg_transcripts_expr), trace = "none", 
          col =  colorRampPalette(c("blue", "white", "red"))(7),
          scale = "row", 
          dendrogram = "row",
          labRow = FALSE,
          Colv = colnames(reg_transcripts_expr),
          Rowv = TRUE,
          hclustfun = function(x) hclust(x, method = "complete"), 
          distfun = function(x) dist(x, method = "euclidean"))
pdf("Results/allgenes_heatmap.pdf", width = 6, height = 6)
print(all_diff_genes_heatmap)
dev.off()
#heatmap 2
reg_transcripts_expr <- avg_vst_data_sub[rownames(avg_vst_data_sub) %in% reg_transcripts$Row,]

#clustering
distanceGene <- dist(reg_transcripts_expr)
distanceSample <- dist(t(reg_transcripts_expr))
clusterGene <- hclust(distanceGene, method="complete")
clusterSample <- hclust(distanceSample, method="average")

reg_transcripts_expr$gene_id <- rownames(reg_transcripts_expr)
deseq2VST <- melt(reg_transcripts_expr, id.vars=c("gene_id"))
deseq2VST$variable <- factor(deseq2VST$variable, levels=clusterSample$labels[clusterSample$order])
reordered_genes <- factor(deseq2VST$gene_id, levels = clusterGene$labels[clusterGene$order])
deseq2VST$gene_id <- reordered_genes
ggplot(deseq2VST, aes(x=variable, y=gene_id, fill=value)) + geom_raster() + scale_fill_viridis(trans="sqrt") + theme(axis.text.x=element_text(angle=65, hjust=1), axis.text.y=element_blank(), axis.ticks.y=element_blank())



#sample distances of averaged data
sampleDists <- dist(t(assay(vst_transformer)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vst_transformer$ID, vst_transformer$sample, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
cluster_heatmap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
pdf("Results/cluster_heatmap.pdf", width = 6, height = 6)
print(cluster_heatmap)
dev.off()

## differential expressed ncRNA only heat map
reg_transcripts_expr_noncoding <- subset(reg_transcripts_expr, rownames(reg_transcripts_expr) %in% noncoding_names)

noncoding_heatmap <- heatmap.2(as.matrix(reg_transcripts_expr_noncoding), trace = "none", 
          col =  colorRampPalette(c("blue", "white", "red"))(7),
          scale = "row", 
          dendrogram = "row",
          labRow = FALSE,
          Colv = colnames(reg_transcripts_expr),
          Rowv = TRUE,
          hclustfun = function(x) hclust(x, method = "complete"), 
          distfun = function(x) dist(x, method = "euclidean"))
pdf("Results/noncoding_heatmap.pdf", width = 6, height = 6)
print(noncoding_heatmap)
dev.off()

## differential expressed coding only heat map
reg_transcripts_expr_coding <- subset(reg_transcripts_expr, rownames(reg_transcripts_expr) %in% coding_names)

coding_heatmap <- heatmap.2(as.matrix(reg_transcripts_expr_coding), trace = "none", 
          col =  colorRampPalette(c("blue", "white", "red"))(7),
          scale = "row", 
          dendrogram = "row",
          labRow = FALSE,
          Colv = colnames(reg_transcripts_expr),
          Rowv = TRUE,
          hclustfun = function(x) hclust(x, method = "complete"), 
          distfun = function(x) dist(x, method = "euclidean"))
pdf("Results/coding_heatmap.pdf", width = 6, height = 6)
print(coding_heatmap)
dev.off()

## calculate stage specificty
calculate_tau <- function(expression_profile) {
  N <- length(expression_profile)
  max_value <- max(expression_profile)
  
  # Normalize expression profile components
  normalized_profile <- expression_profile / max_value
  
  # Calculate τ using the formula
  tau <- sum((1 - normalized_profile)^(N - 1))
  
  return(tau)
}
tau_values <- numeric()
for (i in 1:nrow(avg_vst_data_sub)) {
  gene_expression <- avg_vst_data_sub[i, ]
  tau <- calculate_tau(gene_expression)
  tau_values <- c(tau_values, tau)
}
cat("Gene | Tau\n")
cat("-----------------\n")
for (i in 1:length(tau_values)) {
  cat(i, "   |   ", tau_values[i], "\n")
}



