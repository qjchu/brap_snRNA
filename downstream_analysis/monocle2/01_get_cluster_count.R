library(Seurat)
setwd("bchi_analysis_res")

# CK3D
metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
CK3D_data <- Read10X(data.dir = "../CK3D_1/CK3D_1/outs/filtered_feature_bc_matrix")
CK3D2_data <- Read10X(data.dir = "../CK3D_2/CK3D_2/outs/filtered_feature_bc_matrix")

for (i in 1:20) {
  # i = 1
  sample_cell = rownames(metadata[metadata$harmony_clusters_res1 %in% c(i),])
  sample_cell_CK3D_1 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK3D_1_") { return(sub("^CK3D_1_", "", x))} }))
  CK3D1_data_sample = CK3D_data[,sample_cell_CK3D_1]
  colnames(CK3D1_data_sample) = paste("CK3D_1_", colnames(CK3D1_data_sample), sep = '')

  sample_cell_CK3D_2 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK3D_2_") { return(sub("^CK3D_2_", "", x))} }))
  CK3D2_data_sample = CK3D2_data[,sample_cell_CK3D_2]
  colnames(CK3D2_data_sample) = paste("CK3D_2_", colnames(CK3D2_data_sample), sep = '')

  data_cluster = cbind(CK3D1_data_sample, CK3D2_data_sample)
  expr_matrix <- as(as.matrix(data_cluster), 'sparseMatrix')
  saveRDS(expr_matrix, file = paste('CK3D_cluster', i, '_count.rds', sep = ''))
}

# CK5D
if(TRUE){
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  CK5D_data <- Read10X(data.dir = "../CK5D_1/CK5D_1/outs/filtered_feature_bc_matrix")
  CK5D2_data <- Read10X(data.dir = "../CK5D_2/CK5D_2/outs/filtered_feature_bc_matrix")

  for (i in 1:18) {
    # i = 1
    sample_cell = rownames(metadata[metadata$harmony_clusters_res1 %in% c(i),])
    sample_cell_CK5D_1 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK5D_1_") { return(sub("^CK5D_1_", "", x))} }))
    CK5D1_data_sample = CK5D_data[,sample_cell_CK5D_1]
    colnames(CK5D1_data_sample) = paste("CK5D_1_", colnames(CK5D1_data_sample), sep = '')

    sample_cell_CK5D_2 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK5D_2_") { return(sub("^CK5D_2_", "", x))} }))
    CK5D2_data_sample = CK5D2_data[,sample_cell_CK5D_2]
    colnames(CK5D2_data_sample) = paste("CK5D_2_", colnames(CK5D2_data_sample), sep = '')

    data_cluster = cbind(CK5D1_data_sample, CK5D2_data_sample)
    expr_matrix <- as(as.matrix(data_cluster), 'sparseMatrix')
    saveRDS(expr_matrix, file = paste('CK5D_cluster', i, '_count.rds', sep = ''))
  }

}

# CK7D
if(TRUE){
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  CK7D_data <- Read10X(data.dir = "../CK7D_1/CK7D_1/outs/filtered_feature_bc_matrix")
  CK7D2_data <- Read10X(data.dir = "../CK/CK_bchi/outs/filtered_feature_bc_matrix")

  for (i in 1:25) {
    # i = 1
    sample_cell = rownames(metadata[metadata$harmony_clusters_res1 %in% c(i),])
    sample_cell_CK7D_1 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK7D_1_") { return(sub("^CK7D_1_", "", x))} }))
    CK7D1_data_sample = CK7D_data[,sample_cell_CK7D_1]
    colnames(CK7D1_data_sample) = paste("CK7D_1_", colnames(CK7D1_data_sample), sep = '')

    sample_cell_CK7D_2 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK7D_2_") { return(sub("^CK7D_2_", "", x))} }))
    CK7D2_data_sample = CK7D2_data[,sample_cell_CK7D_2]
    colnames(CK7D2_data_sample) = paste("CK7D_2_", colnames(CK7D2_data_sample), sep = '')

    data_cluster = cbind(CK7D1_data_sample, CK7D2_data_sample)
    expr_matrix <- as(as.matrix(data_cluster), 'sparseMatrix')
    saveRDS(expr_matrix, file = paste('CK7D_cluster', i, '_count.rds', sep = ''))
  }

}

# CK9D
if(TRUE){
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  CK9D_data <- Read10X(data.dir = "../CK9D_1/CK9D_1/outs/filtered_feature_bc_matrix")
  CK9D2_data <- Read10X(data.dir = "../CK9D_2/CK9D_2/outs/filtered_feature_bc_matrix")

  for (i in 1:23) {
    # i = 1
    sample_cell = rownames(metadata[metadata$harmony_clusters_res1 %in% c(i),])
    sample_cell_CK9D_1 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK9D_1_") { return(sub("^CK9D_1_", "", x))} }))
    CK9D1_data_sample = CK9D_data[,sample_cell_CK9D_1]
    colnames(CK9D1_data_sample) = paste("CK9D_1_", colnames(CK9D1_data_sample), sep = '')

    sample_cell_CK9D_2 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 7) == "CK9D_2_") { return(sub("^CK9D_2_", "", x))} }))
    CK9D2_data_sample = CK9D2_data[,sample_cell_CK9D_2]
    colnames(CK9D2_data_sample) = paste("CK9D_2_", colnames(CK9D2_data_sample), sep = '')

    data_cluster = cbind(CK9D1_data_sample, CK9D2_data_sample)
    expr_matrix <- as(as.matrix(data_cluster), 'sparseMatrix')
    saveRDS(expr_matrix, file = paste('CK9D_cluster', i, '_count.rds', sep = ''))
  }

}

# CK11D
if(TRUE){
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  CK11D_data <- Read10X(data.dir = "../CK11D_1/CK11D_1/outs/filtered_feature_bc_matrix")
  CK11D2_data <- Read10X(data.dir = "../CK11D_2/CK11D_2/outs/filtered_feature_bc_matrix")

  for (i in 1:27) {
    # i = 1
    sample_cell = rownames(metadata[metadata$harmony_clusters_res1 %in% c(i),])
    sample_cell_CK11D_1 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 8) == "CK11D_1_") { return(sub("^CK11D_1_", "", x))} }))
    CK11D1_data_sample = CK11D_data[,sample_cell_CK11D_1]
    colnames(CK11D1_data_sample) = paste("CK11D_1_", colnames(CK11D1_data_sample), sep = '')

    sample_cell_CK11D_2 = unlist(sapply(sample_cell, function(x) { if (substr(x, 1, 8) == "CK11D_2_") { return(sub("^CK11D_2_", "", x))} }))
    CK11D2_data_sample = CK11D2_data[,sample_cell_CK11D_2]
    colnames(CK11D2_data_sample) = paste("CK11D_2_", colnames(CK11D2_data_sample), sep = '')

    data_cluster = cbind(CK11D1_data_sample, CK11D2_data_sample)
    expr_matrix <- as(as.matrix(data_cluster), 'sparseMatrix')
    saveRDS(expr_matrix, file = paste('CK11D_cluster', i, '_count.rds', sep = ''))
  }

}
