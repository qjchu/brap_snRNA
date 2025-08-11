library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

sessionInfo()

if(TRUE){
  TT5D_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT5D_1/TT5D_1/outs/filtered_feature_bc_matrix")
  TT5D_1 <- CreateSeuratObject(counts = TT5D_data, project = "TT5D_1", min.cells = 3, min.features = 200)
  
  sce = TT5D_1
  head(sce@meta.data)
  
  q = quantile(sce@meta.data[["nFeature_RNA"]], probs = c(0, 0.05, 0.95, 1))
  sce <- subset(sce, subset = nFeature_RNA > q[2] & nFeature_RNA < q[3])
  
  sce
  
  reticulate::use_python("/public4/chuqj/software/miniconda3/envs/R4/bin/python")
  reticulate::py_module_available("leidenalg")
  
  # run standard anlaysis workflow
  sce <- NormalizeData(sce)
  sce <- FindVariableFeatures(sce)
  sce <- ScaleData(sce)
  sce <- RunPCA(sce)
  sce <- FindNeighbors(sce, dims = 1:50, reduction = "pca")
  sce <- FindClusters(sce, resolution = 1, method="igraph", cluster.name = "unintegrated_clusters_res1", algorithm = 4)
  sce <- RunUMAP(sce, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
  sce <- RunTSNE(sce, dims = 1:50, reduction = "pca", reduction.name = "tsne.unintegrated")
  
  saveRDS(sce, file = "TT5D_merge_res1.rds")
  
  sce@meta.data$umapx.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$umapy.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$tsney.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,2]
  
  write.table(sce@meta.data, file = "TT5D_metadata_res1.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  sce@meta.data$seurat_clusters = sce@meta.data$unintegrated_clusters_res1
  Idents(sce) = sce@meta.data$unintegrated_clusters_res1
  levels(sce)
  
  # cluster markers
  # sce = PrepSCTFindMarkers(sce)
  sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  write.table(sce.markers, file = "TT5D_res1_cluster_markers.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # orig.ident sum
  sample_averages <- AggregateExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_origident_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # unintegrated_clusters_res1 sum
  sample_averages <- AggregateExpression(sce, group.by = "unintegrated_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_clusters_res1_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # orig.ident average
  sample_averages <- AverageExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_origident_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # unintegrated_clusters_res1 average
  sample_averages <- AverageExpression(sce, group.by = "unintegrated_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_clusters_res1_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
}
