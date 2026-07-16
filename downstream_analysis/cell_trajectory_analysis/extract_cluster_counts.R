###############################################################################
# Purpose: Extract raw count matrices for each cluster from a Seurat object
#          and save them as separate RDS files. This is useful for downstream
#          analyses (e.g., pseudobulk DE, gene set enrichment) at the cluster level.
# Input:  - Metadata file containing cluster assignments (harmony_clusters_res1)
#         - Cell Ranger count matrices for each sample (two replicates per timepoint)
# Output: One RDS file per cluster (sparse matrix) named "CK3D_clusterX_count.rds"
###############################################################################

library(Seurat)

# ----------------------------------------------------------------------------
# 1. Set working directory (modify path as needed)
# ----------------------------------------------------------------------------
setwd("path/to/your/analysis_directory")   # e.g., "bchi_analysis_res"

# ----------------------------------------------------------------------------
# 2. Load metadata (contains cluster IDs for each cell)
# ----------------------------------------------------------------------------
metadata <- read.table("CK3D_metadata_res1.txt", header = TRUE)

# ----------------------------------------------------------------------------
# 3. Load 10x count matrices for each replicate (CK3D_1 and CK3D_2)
#    Adjust paths to your Cell Ranger output folders.
# ----------------------------------------------------------------------------
counts_rep1 <- Read10X(data.dir = "path/to/CK3D_1/outs/filtered_feature_bc_matrix")
counts_rep2 <- Read10X(data.dir = "path/to/CK3D_2/outs/filtered_feature_bc_matrix")

# ----------------------------------------------------------------------------
# 4. Loop over clusters (here clusters 1 to 20; adjust as needed)
# ----------------------------------------------------------------------------
for (i in 1:20) {
  
  # 4a. Get cell barcodes for cluster i (from metadata)
  cluster_cells <- rownames(metadata[metadata$harmony_clusters_res1 == i, ])
  
  # 4b. Extract barcodes for each replicate by stripping the sample prefix
  #     Metadata barcodes are like "CK3D_1_AAACGG..."; keep only the suffix.
  barcodes_rep1 <- unlist(sapply(cluster_cells, function(x) {
    if (substr(x, 1, 7) == "CK3D_1_") {
      return(sub("^CK3D_1_", "", x))
    }
  }))
  
  barcodes_rep2 <- unlist(sapply(cluster_cells, function(x) {
    if (substr(x, 1, 7) == "CK3D_2_") {
      return(sub("^CK3D_2_", "", x))
    }
  }))
  
  # 4c. Subset the count matrices for these barcodes
  mat_rep1 <- counts_rep1[, barcodes_rep1, drop = FALSE]
  mat_rep2 <- counts_rep2[, barcodes_rep2, drop = FALSE]
  
  # 4d. Restore full barcodes (with sample prefix) as column names
  colnames(mat_rep1) <- paste0("CK3D_1_", colnames(mat_rep1))
  colnames(mat_rep2) <- paste0("CK3D_2_", colnames(mat_rep2))
  
  # 4e. Combine both replicates into one sparse matrix
  cluster_mat <- cbind(mat_rep1, mat_rep2)
  cluster_sparse <- as(as.matrix(cluster_mat), "sparseMatrix")
  
  # 4f. Save to RDS file
  saveRDS(cluster_sparse, file = paste0("CK3D_cluster", i, "_count.rds"))
}
