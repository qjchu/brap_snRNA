###############################################################################
# Description:
#   - Loads two 10x Genomics snRNA-seq samples (TT7D_1 and TT7D_2)
#   - Merges them into a single Seurat object
#   - Performs quality filtering based on nFeature_RNA
#   - Runs standard Seurat workflow (normalization, PCA, clustering, UMAP/tSNE)
#   - Applies three integration methods: CCA, RPCA, Harmony
#   - Saves the integrated object and exports metadata, marker genes, and aggregated expression
###############################################################################

# Load required libraries
library(dplyr)          # Data manipulation
library(Seurat)         # Single-cell analysis
library(patchwork)      # Plot combination
library(ggplot2)        # Plotting
library(svglite)        # SVG export (not actively used in this script)

# Print session information for reproducibility
sessionInfo()

# Main analysis block (wrapped in if(TRUE) for easy toggling)
if(TRUE){
  # Set working directory to the analysis results folder
  setwd("/public2/chuqj/brap_snRNA/bchi_analysis_res")
  
  # ------------------- Step 1: Load data from Cell Ranger output -------------------
  # Read 10x feature-barcode matrices for sample TT7D_1
  TT7D_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT7D_1/TT7D_1/outs/filtered_feature_bc_matrix")
  # Create a Seurat object for sample 1; keep genes detected in at least 3 cells and cells with at least 200 genes
  TT7D_1 <- CreateSeuratObject(counts = TT7D_data, project = "TT7D_1", min.cells = 3, min.features = 200)
  
  # Read matrix for sample TT7D_2
  TT7D2_data <- Read10X(data.dir = "/public2/chuqj/brap_snRNA/TT7D_2/TT7D_2/outs/filtered_feature_bc_matrix")
  TT7D_2 <- CreateSeuratObject(counts = TT7D2_data, project = "TT7D_2", min.cells = 3, min.features = 200)
  
  # ------------------- Step 2: Merge the two samples -------------------
  # Merge objects; add.cell.ids prefixes each barcode with the sample name
  sce = merge(TT7D_1, y = c(TT7D_2), add.cell.ids = c('TT7D_1', 'TT7D_2'), project = 'merge')
  # View the first few rows of metadata to confirm merging
  head(sce@meta.data)
  
  # ------------------- Step 3: Quality filtering based on number of detected genes -------------------
  # Compute quantiles of nFeature_RNA (0%, 5%, 95%, 100%)
  q = quantile(sce@meta.data[["nFeature_RNA"]], probs = c(0, 0.05, 0.95, 1))
  # Keep cells with nFeature_RNA between the 5th and 95th percentiles (remove outliers)
  sce <- subset(sce, subset = nFeature_RNA > q[2] & nFeature_RNA < q[3])
  
  # The following line is commented because layers are already split; it would split RNA assay by orig.ident
  # sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  # Display object summary
  sce
  
  # ------------------- Step 4: Set Python environment for Leiden algorithm -------------------
  # Specify the Python executable that has leidenalg installed
  reticulate::use_python("/public4/chuqj/software/miniconda3/envs/R4/bin/python")
  # Check that leidenalg is available (will raise error if not)
  reticulate::py_module_available("leidenalg")
  
  # ------------------- Step 5: Standard Seurat workflow (unintegrated) -------------------
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
  
  # ------------------- Step 6: Integration using CCA (Canonical Correlation Analysis) -------------------
  # Integrate layers using CCA; creates a new reduction "integrated.cca"
  sce <- IntegrateLayers(object = sce, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
  # Join the split layers back into a single RNA assay
  sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # Neighbour finding on integrated CCA space
  sce <- FindNeighbors(sce, reduction = "integrated.cca", dims = 1:50)
  # Clustering on CCA space
  sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "cca_clusters_res1")
  # UMAP and t-SNE on CCA reduction
  sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.cca", reduction.name = "umap.cca")
  sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.cca", reduction.name = "tsne.cca")
  
  # ------------------- Step 7: Integration using RPCA (Reciprocal PCA) -------------------
  # Split the RNA assay by sample (required for RPCA integration)
  sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  # Run RPCA integration
  sce <- IntegrateLayers(object = sce, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca", verbose = FALSE)
  # Join layers back
  sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # Neighbour finding, clustering, UMAP, t-SNE on RPCA space
  sce <- FindNeighbors(sce, reduction = "integrated.rpca", dims = 1:50)
  sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "rpca_clusters_res1")
  sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.rpca", reduction.name = "umap.rpca")
  sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.rpca", reduction.name = "tsne.rpca")
  
  # ------------------- Step 8: Integration using Harmony -------------------
  # Split RNA assay again for Harmony (similar requirement)
  sce[["RNA"]] <- split(sce[["RNA"]], f = sce$orig.ident)
  # Run Harmony integration
  sce <- IntegrateLayers(object = sce, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony", verbose = FALSE)
  # Join layers
  sce[["RNA"]] <- JoinLayers(sce[["RNA"]])
  # Neighbour finding, clustering, UMAP, t-SNE on Harmony space
  sce <- FindNeighbors(sce, reduction = "integrated.harmony", dims = 1:50)
  sce <- FindClusters(sce, resolution = 1, method="igraph", algorithm = 4, cluster.name = "harmony_clusters_res1")
  sce <- RunUMAP(sce, dims = 1:50, reduction = "integrated.harmony", reduction.name = "umap.harmony")
  sce <- RunTSNE(sce, dims = 1:50, reduction = "integrated.harmony", reduction.name = "tsne.harmony")
  
  # ------------------- Step 9: Save the integrated Seurat object -------------------
  saveRDS(sce, file = "TT7D_merge_res1.rds")
  
  # ------------------- Step 10: Extract embedding coordinates into metadata for easy export -------------------
  # Unintegrated
  sce@meta.data$umapx.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$umapy.unintegrated = sce@reductions[["umap.unintegrated"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,1]
  sce@meta.data$tsney.unintegrated = sce@reductions[["tsne.unintegrated"]]@cell.embeddings[,2]
  
  # CCA
  sce@meta.data$umapx.cca = sce@reductions[["umap.cca"]]@cell.embeddings[,1]
  sce@meta.data$umapy.cca = sce@reductions[["umap.cca"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.cca = sce@reductions[["tsne.cca"]]@cell.embeddings[,1]
  sce@meta.data$tsney.cca = sce@reductions[["tsne.cca"]]@cell.embeddings[,2]
  
  # RPCA
  sce@meta.data$umapx.rpca = sce@reductions[["umap.rpca"]]@cell.embeddings[,1]
  sce@meta.data$umapy.rpca = sce@reductions[["umap.rpca"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.rpca = sce@reductions[["tsne.rpca"]]@cell.embeddings[,1]
  sce@meta.data$tsney.rpca = sce@reductions[["tsne.rpca"]]@cell.embeddings[,2]
  
  # Harmony
  sce@meta.data$umapx.harmony = sce@reductions[["umap.harmony"]]@cell.embeddings[,1]
  sce@meta.data$umapy.harmony = sce@reductions[["umap.harmony"]]@cell.embeddings[,2]
  sce@meta.data$tsnex.harmony = sce@reductions[["tsne.harmony"]]@cell.embeddings[,1]
  sce@meta.data$tsney.harmony = sce@reductions[["tsne.harmony"]]@cell.embeddings[,2]
  
  # Write the complete metadata table to a text file
  write.table(sce@meta.data, file = "TT7D_metadata_res1.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 11: Set default identity to Harmony clusters (for downstream analysis) -------------------
  sce@meta.data$seurat_clusters = sce@meta.data$harmony_clusters_res1
  Idents(sce) = sce@meta.data$harmony_clusters_res1
  # Show the cluster levels
  levels(sce)
  
  # ------------------- Step 12: Find all marker genes for each cluster (using Harmony clusters) -------------------
  # Note: PrepSCTFindMarkers is commented out because we use RNA assay, not SCT
  # sce = PrepSCTFindMarkers(sce)
  sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  # Write markers to file
  write.table(sce.markers, file = "TT7D_res1_harmony_cluster_markers.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 13: Export aggregated expression (sum) by sample and by cluster -------------------
  # Sum expression per gene across cells from each original sample
  sample_averages <- AggregateExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT7D_RNA_origident_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # Sum expression per gene across cells in each Harmony cluster
  sample_averages <- AggregateExpression(sce, group.by = "harmony_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT7D_RNA_harmony_clusters_res1_sum.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # ------------------- Step 14: Export average expression (mean) by sample and by cluster -------------------
  # Average expression per gene across cells from each sample
  sample_averages <- AverageExpression(sce, group.by = "orig.ident")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT7D_RNA_origident_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
  
  # Average expression per gene across cells in each Harmony cluster
  sample_averages <- AverageExpression(sce, group.by = "harmony_clusters_res1")
  rna_averages <- sample_averages$RNA
  head(rna_averages)
  gene_tpm <- data.matrix(rna_averages)
  head(gene_tpm)
  write.table(gene_tpm, file = "TT7D_RNA_harmony_clusters_res1_avg.txt", sep = '\t', row.names = TRUE, col.names = TRUE)
}
