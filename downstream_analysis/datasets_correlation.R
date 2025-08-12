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

# PNAS LCM corr
if(TRUE){
  exprSet = read.table("GSE12404_expr_matrix_normalized.txt", sep="\t", header = TRUE)
  # transfer to bchi
  if(TRUE){
    transfer = read.table('atha_to_bchi_probe_80.txt', sep = '\t')
    colnames(transfer) = c('probe', 'atha_id', 'bchi_id')
    rownames(transfer) = transfer$probe
    transfer$bchi_id = gsub("\\..*", "", transfer$bchi_id)
    
    exprSet$bchi_id = transfer[rownames(exprSet),]$bchi_id
    exprSet = exprSet[exprSet$bchi_id != '',]
    unique_last_col_indices <- !duplicated(exprSet$bchi_id) & !duplicated(exprSet$bchi_id, fromLast = TRUE)
    exprSet = exprSet[unique_last_col_indices,]
    rownames(exprSet) = exprSet$bchi_id
    exprSet = exprSet[,-88]
    
    # exprSet[, sapply(exprSet, is.character)] <- sapply(exprSet[, sapply(exprSet, is.character)], as.numeric)
  }
  max(exprSet)
  min(exprSet)
  
  exprSet2 <- read.table('TT9D_RNA_harmony_clusters_res1_avg.txt', header = TRUE)
  colnames(exprSet2)
  rownames(exprSet2)
  exprSet2 = exprSet2[rowSums(exprSet2) > 1,]
  
  exprSet2 = exprSet2[rownames(exprSet),]
  exprSet2 = na.omit(exprSet2)
  
  max(exprSet2)
  min(exprSet2)
  
  dat = cbind(exprSet[rownames(exprSet2),], exprSet2)
  colnames(dat)
  
  dat_temp = dat[,c(1:27,30:53,56:67,74:85,88:109)]
  
  # no batch removing
  if(TRUE){
    a = scale(dat_temp)
    
    # cor
    M3 <- cor(a, method = 'spearman') 
    pheatmap::pheatmap(M3)
    
    # clustertree
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    # pca
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    class(dat_pca)
    fviz_pca_ind(dat_pca)
    
    # anthor pca
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
  
  # remove batch -- combat
  if(TRUE){
    a = ComBat(scale(dat_temp), batch = c(rep('batch1',75), rep('batch2',22)))
    
    # cor
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3[76:97,1:75])
    
    # clustertree
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    # pca
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    class(dat_pca)
    fviz_pca_ind(dat_pca)
    
    # anthor pca
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
}

# NaturePlant Corr
if(TRUE){
  exprSet <- read.table('GSE157145_CPM_average_expression.txt', header = TRUE)
  colnames(exprSet)
  rownames(exprSet) = exprSet[,1]
  exprSet = exprSet[,-1]
  
  # transfer to bchi
  if(TRUE){
    transfer = read.table('atha_to_bchi_besthit_80.txt', sep = ' ')
    colnames(transfer) = c('atha_id', 'bchi_id')
    # transfer = unique(transfer)
    transfer = transfer[!duplicated(transfer$atha_id), ]
    rownames(transfer) = transfer$atha_id
    transfer$bchi_id = gsub("\\..*", "", transfer$bchi_id)
    
    exprSet$bchi_id = transfer[rownames(exprSet),]$bchi_id
    exprSet = exprSet[exprSet$bchi_id != '',]
    unique_last_col_indices <- !duplicated(exprSet$bchi_id) & !duplicated(exprSet$bchi_id, fromLast = TRUE)
    exprSet = exprSet[unique_last_col_indices,]
    rownames(exprSet) = exprSet$bchi_id
    exprSet = exprSet[,-40]
    
    # exprSet[, sapply(exprSet, is.character)] <- sapply(exprSet[, sapply(exprSet, is.character)], as.numeric)
  }
  max(exprSet)
  min(exprSet)
  
  M <- cor(exprSet, method = 'spearman')
  pheatmap::pheatmap(M)
  
  exprSet2 <- read.table('../TT9D/TT9D_RNA_harmony_clusters_res1_avg.txt', header = TRUE)
  colnames(exprSet2)
  rownames(exprSet2)
  exprSet2 = exprSet2[rowSums(exprSet2) > 1,]
  
  exprSet2 = exprSet2[rownames(exprSet),]
  exprSet2 = na.omit(exprSet2)
  
  max(exprSet2)
  min(exprSet2)
  
  dat = cbind(exprSet[rownames(exprSet2),], exprSet2)
  colnames(dat)
  
  dat_temp = dat
  
  # no batch removing
  if(TRUE){
    a = scale(dat_temp)
    
    # cor
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3)
    
    # clustertree
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    # pca
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    class(dat_pca)
    fviz_pca_ind(dat_pca)
    
    # anthor pca
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
  
  # remove batch -- combat
  if(TRUE){
    a = ComBat(scale(dat_temp), batch = c(rep('batch1',39), rep('batch2',27)))
    
    # cor
    M3 <- cor(a, method = 'spearman')
    pheatmap::pheatmap(M3[40:66,1:39])
    pheatmap::pheatmap(M3[40:66,1:21])
    pheatmap::pheatmap(M3[40:66,22:39])
    
    # clustertree
    sampleTree <- hclust(dist(t(a)), method = "average")
    plot(sampleTree)
    
    # pca
    b <- as.data.frame(t(a))
    dat_pca <- PCA(b[,-ncol(b)], graph = FALSE)
    class(dat_pca)
    fviz_pca_ind(dat_pca)
    
    # anthor pca
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
}
