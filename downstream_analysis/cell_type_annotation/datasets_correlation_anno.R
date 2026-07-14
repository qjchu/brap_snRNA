###############################################################################
# Script: Cross-dataset correlation analysis between external LCM/RNA-seq data
#         and our snRNA-seq cluster averages (example: TT9D).
# Steps for each dataset:
#   1. Load external expression data and convert Arabidopsis gene IDs to B. rapa IDs.
#   2. Load our cluster-averaged expression matrix.
#   3. Merge both matrices (common genes).
#   4. Scale data and perform correlation, clustering, PCA before/after ComBat correction.
###############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(sva)
library(ggsci)
library(ggrepel)
library(tidyr)
library(stringr)

sessionInfo()
getwd()

# ------------------- Section 1: PNAS LCM dataset correlation -------------------
if(TRUE){
  # 1a. Read external LCM expression matrix (GSE12404)
  exprSet = read.table("GSE12404_expr_matrix_normalized.txt", sep="\t", header = TRUE)
  
  # 1b. Convert Arabidopsis probe IDs to B. rapa gene IDs using a mapping file
  if(TRUE){
    transfer = read.table('atha_to_bchi_probe_80.txt', sep = '\t')
    colnames(transfer) = c('probe', 'atha_id', 'bchi_id')
    rownames(transfer) = transfer$probe
    transfer$bchi_id = gsub("\\..*", "", transfer$bchi_id)   # Strip version suffix
    
    exprSet$bchi_id = transfer[rownames(exprSet),]$bchi_id
    exprSet = exprSet[exprSet$bchi_id != '',]                # Remove unmapped probes
    # Keep only genes with unique mapping (no duplicates)
    unique_last_col_indices <- !duplicated(exprSet$bchi_id) & !duplicated(exprSet$bchi_id, fromLast = TRUE)
    exprSet = exprSet[unique_last_col_indices,]
    rownames(exprSet) = exprSet$bchi_id
    exprSet = exprSet[,-88]                                   # Remove the ID column
  }
  max(exprSet)
  min(exprSet)
  
  # 1c. Read our cluster-averaged expression (TT9D)
  exprSet2 <- read.table('TT9D_RNA_harmony_clusters_res1_avg.txt', header = TRUE)
  colnames(exprSet2)
  rownames(exprSet2)
  exprSet2 = exprSet2[rowSums(exprSet2) > 1, ]               # Remove genes with low total expression
  
  # Keep only overlapping genes
  exprSet2 = exprSet2[rownames(exprSet), ]
  exprSet2 = na.omit(exprSet2)
  
  max(exprSet2)
  min(exprSet2)
  
  # Merge external and our data
  dat = cbind(exprSet[rownames(exprSet2), ], exprSet2)
  colnames(dat)
  
  # Select specific columns (adjust indices to keep relevant samples)
  dat_temp = dat[, c(1:27, 30:53, 56:67, 74:85, 88:109)]
  
  # 1d. Analysis without batch correction
  if(TRUE){
    a = scale(dat_temp)   # Scale genes (columns) to zero mean and unit variance
    
    # Correlation heatmap
    M3 <- cor(a, method = 'spearman') 
    pheatmap::pheatmap(M3)
    
    # Hierarchical clustering of samples
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    # PCA using factoextra (plots not shown here)
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    class(dat_pca)
    fviz_pca_ind(dat_pca)
    
    # Alternative PCA with prcomp and labelled plot
    pca_dat <- prcomp(a, scale. = FALSE)
    head(pca_dat$rotation)
    pca.data <- data.frame(sample = rownames(pca_dat$rotation),
                           Type = sapply(rownames(pca_dat$rotation), function(x) gsub("\\..*", "", x)), 
                           pca_dat$rotation)
    ggscatter(pca.data, x = "PC1", y = "PC2", color = "Type", label = 'sample', repel = TRUE) + 
      theme(text =  element_text(size = 12),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            legend.position = "right")
  }
  
  # 1e. Analysis with ComBat batch correction (removes batch effects between datasets)
  if(TRUE){
    # Batch vector: first 75 columns from external data, last 22 from ours
    a = ComBat(scale(dat_temp), batch = c(rep('batch1', 75), rep('batch2', 22)))
    
    # Repeat correlation, clustering, PCA on corrected data
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3[76:97, 1:75])   # Cross-correlation between batches
    
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    fviz_pca_ind(dat_pca)
    
    pca_dat <- prcomp(a, scale. = FALSE)
    pca.data <- data.frame(sample = rownames(pca_dat$rotation),
                           Type = sapply(rownames(pca_dat$rotation), function(x) gsub("\\..*", "", x)), 
                           pca_dat$rotation)
    ggscatter(pca.data, x = "PC1", y = "PC2", color = "Type", label = 'sample', repel = TRUE) + 
      theme(text =  element_text(size = 12),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            legend.position = "right")
  }
}

# ------------------- Section 2: NaturePlant RNA-seq dataset correlation -------------------
if(TRUE){
  # 2a. Read external CPM average expression (GSE157145)
  exprSet <- read.table('GSE157145_CPM_average_expression.txt', header = TRUE)
  colnames(exprSet)
  rownames(exprSet) = exprSet[,1]
  exprSet = exprSet[,-1]
  
  # 2b. Convert Arabidopsis gene IDs to B. rapa (using a different mapping file)
  if(TRUE){
    transfer = read.table('atha_to_bchi_besthit_80.txt', sep = ' ')
    colnames(transfer) = c('atha_id', 'bchi_id')
    transfer = transfer[!duplicated(transfer$atha_id), ]    # Keep only unique Arabidopsis IDs
    rownames(transfer) = transfer$atha_id
    transfer$bchi_id = gsub("\\..*", "", transfer$bchi_id)
    
    exprSet$bchi_id = transfer[rownames(exprSet), ]$bchi_id
    exprSet = exprSet[exprSet$bchi_id != '', ]
    unique_last_col_indices <- !duplicated(exprSet$bchi_id) & !duplicated(exprSet$bchi_id, fromLast = TRUE)
    exprSet = exprSet[unique_last_col_indices, ]
    rownames(exprSet) = exprSet$bchi_id
    exprSet = exprSet[,-40]   # Remove ID column (40th column)
  }
  max(exprSet)
  min(exprSet)
  
  # Correlation within external dataset
  M <- cor(exprSet, method = 'spearman')
  pheatmap::pheatmap(M)
  
  # 2c. Load our TT9D average expression (same as before)
  exprSet2 <- read.table('../TT9D/TT9D_RNA_harmony_clusters_res1_avg.txt', header = TRUE)
  colnames(exprSet2)
  rownames(exprSet2)
  exprSet2 = exprSet2[rowSums(exprSet2) > 1, ]
  
  exprSet2 = exprSet2[rownames(exprSet), ]
  exprSet2 = na.omit(exprSet2)
  
  max(exprSet2)
  min(exprSet2)
  
  # Merge datasets
  dat = cbind(exprSet[rownames(exprSet2), ], exprSet2)
  colnames(dat)
  
  dat_temp = dat   # All columns used (no subsetting here)
  
  # 2d. Without batch correction
  if(TRUE){
    a = scale(dat_temp)
    
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3)
    
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    fviz_pca_ind(dat_pca)
    
    pca_dat <- prcomp(a, scale. = FALSE)
    pca.data <- data.frame(sample = rownames(pca_dat$rotation),
                           Type = sapply(rownames(pca_dat$rotation), function(x) gsub("\\..*", "", x)), 
                           pca_dat$rotation)
    ggscatter(pca.data, x = "PC1", y = "PC2", color = "Type", label = 'sample', repel = TRUE) + 
      theme(text =  element_text(size = 12),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            legend.position = "right")
  }
  
  # 2e. With ComBat batch correction (39 external samples vs 27 clusters)
  if(TRUE){
    a = ComBat(scale(dat_temp), batch = c(rep('batch1', 39), rep('batch2', 27)))
    
    # Explore cross-correlation matrices for subsets
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3[40:66, 1:39])
    pheatmap::pheatmap(M3[40:66, 1:21])
    pheatmap::pheatmap(M3[40:66, 22:39])
    
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    fviz_pca_ind(dat_pca)
    
    pca_dat <- prcomp(a, scale. = FALSE)
    pca.data <- data.frame(sample = rownames(pca_dat$rotation),
                           Type = sapply(rownames(pca_dat$rotation), function(x) gsub("\\..*", "", x)), 
                           pca_dat$rotation)
    ggscatter(pca.data, x = "PC1", y = "PC2", color = "Type", label = 'sample', repel = TRUE) + 
      theme(text =  element_text(size = 12),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            legend.position = "right")
  }
}
