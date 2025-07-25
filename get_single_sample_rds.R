library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

sessionInfo()

if(TRUE){
  TT5D_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT5D_1/TT5D_1/outs/filtered_feature_bc_matrix")
  TT5D_1 <- CreateSeuratObject(counts = TT5D_data, project = "TT5D_1", min.cells = 3, min.features = 200)
  
  # TT5D2_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT5D_2/TT5D_2/outs/filtered_feature_bc_matrix")
  # TT5D_2 <- CreateSeuratObject(counts = TT5D2_data, project = "TT5D_2", min.cells = 3, min.features = 200)
  
  # sce = merge(TT5D_1, y = c(TT5D_2), add.cell.ids = c('TT5D_1', 'TT5D_2'), project = 'merge')
  sce = TT5D_1
  head(sce@meta.data)
  
  q = quantile(sce@meta.data[["nFeature_RNA"]], probs = c(0, 0.05, 0.95, 1))
  sce <- subset(sce, subset = nFeature_RNA > q[2] & nFeature_RNA < q[3])
  
  # sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident) # The following layers are already split
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
  
  # # cca
  # sce <- IntegrateLayers(object = sce, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
  # sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # sce <- FindNeighbors(sce, reduction = "integrated.cca", dims = 1:50)
  # sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "cca_clusters_res1")
  # sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.cca", reduction.name = "umap.cca")
  # sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.cca", reduction.name = "tsne.cca")
  # 
  # # rpca
  # # no applicable method for 'Assays' applied to an object of class "NULL"
  # sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  # sce <- IntegrateLayers(object = sce, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE)
  # sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # sce <- FindNeighbors(sce, reduction = "integrated.rpca", dims = 1:50)
  # sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "rpca_clusters_res1")
  # sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.rpca", reduction.name = "umap.rpca")
  # sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.rpca", reduction.name = "tsne.rpca")
  # 
  # # harmony
  # # Error in names(groups) <- "group" : attempt to set an attribute on NULL
  # sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  # sce <- IntegrateLayers(object = sce, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony", verbose = FALSE)
  # sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # sce <- FindNeighbors(sce, reduction = "integrated.harmony", dims = 1:50)
  # sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "harmony_clusters_res1")
  # sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.harmony", reduction.name = "umap.harmony")
  # sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.harmony", reduction.name = "tsne.harmony")
  
  saveRDS(sce, file = "TT5D_merge_res1.rds")
  
  sce@meta.data$umapx.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$umapy.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$tsney.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,2]
  
  # sce@meta.data$umapx.cca = sce@reductions[["umap.cca"]]@cell.embeddings[,1]
  # sce@meta.data$umapy.cca = sce@reductions[["umap.cca"]]@cell.embeddings[,2]
  # sce@meta.data$tsnex.cca = sce@reductions[["tsne.cca"]]@cell.embeddings[,1]
  # sce@meta.data$tsney.cca = sce@reductions[["tsne.cca"]]@cell.embeddings[,2]
  # 
  # sce@meta.data$umapx.rpca = sce@reductions[["umap.rpca"]]@cell.embeddings[,1]
  # sce@meta.data$umapy.rpca = sce@reductions[["umap.rpca"]]@cell.embeddings[,2]
  # sce@meta.data$tsnex.rpca = sce@reductions[["tsne.rpca"]]@cell.embeddings[,1]
  # sce@meta.data$tsney.rpca = sce@reductions[["tsne.rpca"]]@cell.embeddings[,2]
  # 
  # sce@meta.data$umapx.harmony = sce@reductions[["umap.harmony"]]@cell.embeddings[,1]
  # sce@meta.data$umapy.harmony = sce@reductions[["umap.harmony"]]@cell.embeddings[,2]
  # sce@meta.data$tsnex.harmony = sce@reductions[["tsne.harmony"]]@cell.embeddings[,1]
  # sce@meta.data$tsney.harmony = sce@reductions[["tsne.harmony"]]@cell.embeddings[,2]
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
