###############################################################################
# Description:
#   - Loads one 10x Genomics snRNA-seq sample (TT5D_1)
#   - Creates a Seurat object and filters cells based on nFeature_RNA
#   - Runs the standard Seurat workflow: normalisation, PCA, clustering, UMAP/tSNE
#   - Saves the processed object, exports metadata, marker genes, and aggregated expression
#   - No integration is performed (single sample)
###############################################################################

# Load required libraries
library(dplyr)          # Data manipulation
library(Seurat)         # Single-cell analysis
library(patchwork)      # Plot combination (not actively used)
library(ggplot2)        # Plotting (not actively used)
library(svglite)        # SVG export (not actively used)

# Print session information for reproducibility
sessionInfo()

# Main analysis block (wrapped in if(TRUE) for easy toggling)
if(TRUE){
  # ------------------- Step 1: Load data from Cell Ranger output -------------------
  # Read the 10x feature-barcode matrix for sample TT5D_1
  TT5D_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT5D_1/TT5D_1/outs/filtered_feature_bc_matrix")
  # Create a Seurat object; keep genes detected in at least 3 cells and cells with at least 200 genes
  TT5D_1 <- CreateSeuratObject(counts = TT5D_data, project = "TT5D_1", min.cells = 3, min.features = 200)
  
  # Assign the object to a variable named 'sce' for consistency
  sce = TT5D_1
  # View the first few rows of metadata
  head(sce@meta.data)
  
  # ------------------- Step 2: Quality filtering based on number of detected genes -------------------
  # Compute quantiles of nFeature_RNA (0%, 5%, 95%, 100%)
  q = quantile(sce@meta.data[["nFeature_RNA"]], probs = c(0, 0.05, 0.95, 1))
  # Keep cells with nFeature_RNA between the 5th and 95th percentiles (remove outliers)
  sce <- subset(sce, subset = nFeature_RNA > q[2] & nFeature_RNA < q[3])
  
  # Display object summary
  sce
  
  # ------------------- Step 3: Set Python environment for Leiden algorithm -------------------
  # Specify the Python executable that has leidenalg installed
  reticulate::use_python("/public4/chuqj/software/miniconda3/envs/R4/bin/python")
  # Check that leidenalg is available (will raise error if not)
  reticulate::py_module_available("leidenalg")
  
  # ------------------- Step 4: Standard Seurat workflow (unintegrated) -------------------
  # Normalise the data using LogNormalize (default)
  sce <- NormalizeData(sce)
  # Identify highly variable features for PCA
  sce <- FindVariableFeatures(sce)
  # Scale the data (regress out unwanted variation if needed; here default)
  sce <- ScaleData(sce)
  # Perform PCA on the scaled data
  sce <- RunPCA(sce)
  # Build a nearest-neighbour graph using the first 50 PCs
  sce <- FindNeighbors(sce, dims = 1:50, reduction = "pca")
  # Cluster the cells using the Leiden algorithm (algorithm = 4) at resolution 1.0
  sce <- FindClusters(sce, resolution = 1, method="igraph", cluster.name = "unintegrated_clusters_res1", algorithm = 4)
  # Compute UMAP embedding from PCA
  sce <- RunUMAP(sce, dims = 1:50, reduction = "pca", reduction.name = "umap.unintegrated")
  # Compute t-SNE embedding from PCA
  sce <- RunTSNE(sce, dims = 1:50, reduction = "pca", reduction.name = "tsne.unintegrated")
  
  # ------------------- Step 5: Save the processed Seurat object -------------------
  saveRDS(sce, file = "TT5D_merge_res1.rds")
  
  # ------------------- Step 6: Extract embedding coordinates into metadata for easy export -------------------
  # UMAP coordinates
  sce@meta.data$umapx.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$umapy.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,2]
  # t-SNE coordinates
  sce@meta.data$tsnex.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$tsney.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,2]
  
  # Write the complete metadata table to a text file
  write.table(sce@meta.data, file = "TT5D_metadata_res1.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 7: Set default identity to the unintegrated clusters -------------------
  sce@meta.data$seurat_clusters = sce@meta.data$unintegrated_clusters_res1
  Idents(sce) = sce@meta.data$unintegrated_clusters_res1
  # Show the cluster levels
  levels(sce)
  
  # ------------------- Step 8: Find all marker genes for each cluster -------------------
  # Note: PrepSCTFindMarkers is commented out because we use the RNA assay, not SCT
  # sce = PrepSCTFindMarkers(sce)
  sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # Write markers to file
  write.table(sce.markers, file = "TT5D_res1_cluster_markers.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 9: Export aggregated expression (sum) by sample and by cluster -------------------
  # Sum expression per gene across all cells (only one sample, but grouped by orig.ident)
  sample_averages <- AggregateExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_origident_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # Sum expression per gene across cells in each cluster (unintegrated)
  sample_averages <- AggregateExpression(sce, group.by = "unintegrated_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_clusters_res1_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 10: Export average expression (mean) by sample and by cluster -------------------
  # Average expression per gene across all cells (grouped by orig.ident)
  sample_averages <- AverageExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_origident_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # Average expression per gene across cells in each cluster
  sample_averages <- AverageExpression(sce, group.by = "unintegrated_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT5D_RNA_clusters_res1_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
}
