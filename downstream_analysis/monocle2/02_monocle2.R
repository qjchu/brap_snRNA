library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

setwd("bchi_analysis_res")

## CK3D
if(TRUE){
  if(TRUE){
    metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
    
    data10 <- readRDS('CK3D_cluster10_count.rds')
    pd10 = metadata[colnames(data10),]
    pd10 = as.data.frame(pd10)
    pd10$celltype = 'SC_Endosperm_c10'
    
    data13 <- readRDS('CK3D_cluster13_count.rds')
    pd13 = metadata[colnames(data13),]
    pd13 = as.data.frame(pd13)
    pd13$celltype = 'SC_c13'
    
    data14 <- readRDS('CK3D_cluster14_count.rds')
    pd14 = metadata[colnames(data14),]
    pd14 = as.data.frame(pd14)
    pd14$celltype = 'SC_Endosperm_c14'
    
    data18 <- readRDS('CK3D_cluster18_count.rds')
    pd18 = metadata[colnames(data18),]
    pd18 = as.data.frame(pd18)
    pd18$celltype = 'PEN_c18'
    
    data20 <- readRDS('CK3D_cluster20_count.rds')
    pd20 = metadata[colnames(data20),]
    pd20 = as.data.frame(pd20)
    pd20$celltype = 'CZE_c20'
    
  }
  
  fd <- as.data.frame(rownames(data14))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd18,pd20)
  
  data = cbind(data18,data20)
  rownames(data) = rownames(data14)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK3D_c18c20_monocle.rds')
  write.table(deg, file = "CK3D_c18c20_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  saveRDS(cds, 'CK3D_c10c13c14_monocle.rds')
  write.table(deg, file = "CK3D_c10c13c14_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  saveRDS(cds, 'CK3D_c10c14c18_monocle.rds')
  write.table(deg, file = "CK3D_c10c14c18_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  saveRDS(cds, 'CK3D_c13c14c18_monocle.rds')
  write.table(deg, file = "CK3D_c13c14c18_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

## CK5D
if(TRUE){
  if(TRUE){
    metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
    
    data3 <- readRDS('CK5D_cluster3_count.rds')
    pd3 = metadata[colnames(data3),]
    pd3 = as.data.frame(pd3)
    pd3$celltype = 'PEN_c3'
    
    data6 <- readRDS('CK5D_cluster6_count.rds')
    pd6 = metadata[colnames(data6),]
    pd6 = as.data.frame(pd6)
    pd6$celltype = 'CZSC_c6'
    
    data14 <- readRDS('CK5D_cluster14_count.rds')
    pd14 = metadata[colnames(data14),]
    pd14 = as.data.frame(pd14)
    pd14$celltype = 'CZE_c14'
    
    data16 <- readRDS('CK5D_cluster16_count.rds')
    pd16 = metadata[colnames(data16),]
    pd16 = as.data.frame(pd16)
    pd16$celltype = 'MCE_c16'
    
    data17 <- readRDS('CK5D_cluster17_count.rds')
    pd17 = metadata[colnames(data17),]
    pd17 = as.data.frame(pd17)
    pd17$celltype = 'CZSC_c17'
    
  }
  
  fd <- as.data.frame(rownames(data14))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd3,pd14,pd16)
  
  data = cbind(data3,data14,data16)
  rownames(data) = rownames(data14)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK5D_c3c14c16_monocle.rds')
  write.table(deg, file = "CK5D_c3c14c16_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  if(TRUE){
    fd <- as.data.frame(rownames(data14))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd6,pd17)
    
    data = cbind(data6,data17)
    rownames(data) = rownames(data14)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK5D_c6c17_monocle.rds')
    write.table(deg, file = "CK3D_c6c17_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
  }
  
}

## CK7D
if(TRUE){
  if(TRUE){
    metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
    
    data3 <- readRDS('CK7D_cluster3_count.rds')
    pd3 = metadata[colnames(data3),]
    pd3 = as.data.frame(pd3)
    pd3$celltype = 'PEN_c3'
    
    data4 <- readRDS('CK7D_cluster4_count.rds')
    pd4 = metadata[colnames(data4),]
    pd4 = as.data.frame(pd4)
    pd4$celltype = 'SC_Endosperm_c4'
    
    data8 <- readRDS('CK7D_cluster8_count.rds')
    pd8 = metadata[colnames(data8),]
    pd8 = as.data.frame(pd8)
    pd8$celltype = 'CZSC_c8'
    
    data9 <- readRDS('CK7D_cluster9_count.rds')
    pd9 = metadata[colnames(data9),]
    pd9 = as.data.frame(pd9)
    pd9$celltype = 'CZSC_c9'
    
    data10 <- readRDS('CK7D_cluster10_count.rds')
    pd10 = metadata[colnames(data10),]
    pd10 = as.data.frame(pd10)
    pd10$celltype = 'c10'
    
    data11 <- readRDS('CK7D_cluster11_count.rds')
    pd11 = metadata[colnames(data11),]
    pd11 = as.data.frame(pd11)
    pd11$celltype = 'EP_Diving_Cell_c11'
    
    data13 <- readRDS('CK7D_cluster13_count.rds')
    pd13 = metadata[colnames(data13),]
    pd13 = as.data.frame(pd13)
    pd13$celltype = 'c13'
    
    data14 <- readRDS('CK7D_cluster14_count.rds')
    pd14 = metadata[colnames(data14),]
    pd14 = as.data.frame(pd14)
    pd14$celltype = 'c14'
    
    data15 <- readRDS('CK7D_cluster15_count.rds')
    pd15 = metadata[colnames(data15),]
    pd15 = as.data.frame(pd15)
    pd15$celltype = 'SC_Endosperm_c15'
    
    data16 <- readRDS('CK7D_cluster16_count.rds')
    pd16 = metadata[colnames(data16),]
    pd16 = as.data.frame(pd16)
    pd16$celltype = 'CZE_c16'
    
    data17 <- readRDS('CK7D_cluster17_count.rds')
    pd17 = metadata[colnames(data17),]
    pd17 = as.data.frame(pd17)
    pd17$celltype = 'EP_Diving_Cell_c17'
    
    data18 <- readRDS('CK7D_cluster18_count.rds')
    pd18 = metadata[colnames(data18),]
    pd18 = as.data.frame(pd18)
    pd18$celltype = 'c18'
    
    data19 <- readRDS('CK7D_cluster19_count.rds')
    pd19 = metadata[colnames(data19),]
    pd19 = as.data.frame(pd19)
    pd19$celltype = 'MCE_c19'
    
    data20 <- readRDS('CK7D_cluster20_count.rds')
    pd20 = metadata[colnames(data20),]
    pd20 = as.data.frame(pd20)
    pd20$celltype = 'c20'
    
    data24 <- readRDS('CK7D_cluster24_count.rds')
    pd24 = metadata[colnames(data24),]
    pd24 = as.data.frame(pd24)
    pd24$celltype = 'CZSC_c24'
    
    data25 <- readRDS('CK7D_cluster25_count.rds')
    pd25 = metadata[colnames(data25),]
    pd25 = as.data.frame(pd25)
    pd25$celltype = 'EP_c25'
    
  }
  
  fd <- as.data.frame(rownames(data11))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd8,pd9,pd24)
  
  data = cbind(data8,data9,data24)
  rownames(data) = rownames(data11)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK7D_c8c9c24_monocle.rds')
  write.table(deg, file = "CK7D_c8c9c24_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  saveRDS(cds, 'CK7D_c3c4c15c16c19_monocle.rds')
  write.table(deg, file = "CK7D_c3c4c15c16c19_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  saveRDS(cds, 'CK7D_c11c17c25_monocle.rds')
  write.table(deg, file = "CK7D_c11c17c25_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

getPD = function(data, celltype){
  pd = metadata[colnames(data),]
  pd = as.data.frame(pd)
  pd$celltype = celltype
  return(pd)
}

## CK9D
if(TRUE){
  if(TRUE){
    metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
    
    data1 <- readRDS('CK9D_cluster1_count.rds')
    pd1 = getPD(data = data1, celltype = 'c1')
    
    data4 <- readRDS('CK9D_cluster4_count.rds')
    pd4 = getPD(data = data4, celltype = 'CZSC_c4')
    
    data9 <- readRDS('CK9D_cluster9_count.rds')
    pd9 = getPD(data = data9, celltype = 'c9')
    
    data10 <- readRDS('CK9D_cluster10_count.rds')
    pd10 = getPD(data = data10, celltype = 'EP_Dividing_Cell_c10')
    
    data11 <- readRDS('CK9D_cluster11_count.rds')
    pd11 = getPD(data = data11, celltype = 'c11')
    
    data14 <- readRDS('CK9D_cluster14_count.rds')
    pd14 = getPD(data = data14, celltype = 'CZSC_c14')
    
    data15 <- readRDS('CK9D_cluster15_count.rds')
    pd15 = getPD(data = data15, celltype = 'c15')
    
    data16 <- readRDS('CK9D_cluster16_count.rds')
    pd16 = getPD(data = data16, celltype = 'c16')
    
    data19 <- readRDS('CK9D_cluster19_count.rds')
    pd19 = getPD(data = data19, celltype = 'c19')
    
    data20 <- readRDS('CK9D_cluster20_count.rds')
    pd20 = getPD(data = data20, celltype = 'c20')
    
    data21 <- readRDS('CK9D_cluster21_count.rds')
    pd21 = getPD(data = data21, celltype = 'EP_c21')
    
    data22 <- readRDS('CK9D_cluster22_count.rds')
    pd22 = getPD(data = data22, celltype = 'CZSC_c22')
  }
  
  fd <- as.data.frame(rownames(data11))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd1,pd9,pd11,pd15,pd16,pd19,pd20)
  
  data = cbind(data1,data9,data11,data15,data16,data19,data20)
  rownames(data) = rownames(data11)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK9D_c1c9c11c15c16c19c20_monocle.rds')
  write.table(deg, file = "CK9D_c1c9c11c15c16c19c20_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  if(TRUE){
    fd <- as.data.frame(rownames(data11))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd10,pd21)
    
    data = cbind(data10,data21)
    rownames(data) = rownames(data11)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK9D_c10c21_monocle.rds')
    write.table(deg, file = "CK9D_c10c21_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  
  if(TRUE){
    fd <- as.data.frame(rownames(data11))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd4,pd14,pd22)
    
    data = cbind(data4,data14,data22)
    rownames(data) = rownames(data11)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK9D_c4c14c22_monocle.rds')
    write.table(deg, file = "CK9D_c4c14c22_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}

## CK11D
if(TRUE){
  if(TRUE){
    metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
    
    data1 <- readRDS('CK11D_cluster1_count.rds')
    data2 <- readRDS('CK11D_cluster2_count.rds')
    data3 <- readRDS('CK11D_cluster3_count.rds')
    data4 <- readRDS('CK11D_cluster4_count.rds')
    data5 <- readRDS('CK11D_cluster5_count.rds')
    data6 <- readRDS('CK11D_cluster6_count.rds')
    data7 <- readRDS('CK11D_cluster7_count.rds')
    data8 <- readRDS('CK11D_cluster8_count.rds')
    data9 <- readRDS('CK11D_cluster9_count.rds')
    data10 <- readRDS('CK11D_cluster10_count.rds')
    data11 <- readRDS('CK11D_cluster11_count.rds')
    data12 <- readRDS('CK11D_cluster12_count.rds')
    data13 <- readRDS('CK11D_cluster13_count.rds')
    data14 <- readRDS('CK11D_cluster14_count.rds')
    data15 <- readRDS('CK11D_cluster15_count.rds')
    data16 <- readRDS('CK11D_cluster16_count.rds')
    data17 <- readRDS('CK11D_cluster17_count.rds')
    data18 <- readRDS('CK11D_cluster18_count.rds')
    data19 <- readRDS('CK11D_cluster19_count.rds')
    data20 <- readRDS('CK11D_cluster20_count.rds')
    data21 <- readRDS('CK11D_cluster21_count.rds')
    data22 <- readRDS('CK11D_cluster22_count.rds')
    data23 <- readRDS('CK11D_cluster23_count.rds')
    data24 <- readRDS('CK11D_cluster24_count.rds')
    data25 <- readRDS('CK11D_cluster25_count.rds')
    data26 <- readRDS('CK11D_cluster26_count.rds')
    data27 <- readRDS('CK11D_cluster27_count.rds')
    
    pd1 = getPD(data = data1, celltype = 'c1')
    pd2 = getPD(data = data2, celltype = 'c2')
    pd3 = getPD(data = data3, celltype = 'c3')
    pd4 = getPD(data = data4, celltype = 'c4')
    pd5 = getPD(data = data5, celltype = 'c5')
    pd6 = getPD(data = data6, celltype = 'c6')
    pd7 = getPD(data = data7, celltype = 'c7')
    pd8 = getPD(data = data8, celltype = 'c8')
    pd9 = getPD(data = data9, celltype = 'CZE_c9')
    pd10 = getPD(data = data10, celltype = 'c10')
    pd11 = getPD(data = data11, celltype = 'c11')
    pd12 = getPD(data = data12, celltype = 'CZE_c12')
    pd13 = getPD(data = data13, celltype = 'c13')
    pd14 = getPD(data = data14, celltype = 'c14')
    pd15 = getPD(data = data15, celltype = 'c15')
    pd16 = getPD(data = data16, celltype = 'c16')
    pd17 = getPD(data = data17, celltype = 'c17')
    pd18 = getPD(data = data18, celltype = 'c18')
    pd19 = getPD(data = data19, celltype = 'c19')
    pd20 = getPD(data = data20, celltype = 'c20')
    pd21 = getPD(data = data21, celltype = 'c21')
    pd22 = getPD(data = data22, celltype = 'c22')
    pd23 = getPD(data = data23, celltype = 'CZE_c23')
    pd24 = getPD(data = data24, celltype = 'c24')
    pd25 = getPD(data = data25, celltype = 'c25')
    pd26 = getPD(data = data26, celltype = 'c26')
    pd27 = getPD(data = data27, celltype = 'c27')
  }
  
  fd <- as.data.frame(rownames(data11))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName

  pd = rbind(pd1,pd2,pd6,pd16,pd19,pd22)
  
  data = cbind(data1,data2,data6,data16,data19,data22)
  rownames(data) = rownames(data11)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK11D_c1c2c6c16c19c22_monocle.rds')
  write.table(deg, file = "CK11D_c1c2c6c16c19c22_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  if(TRUE){
    fd <- as.data.frame(rownames(data11))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd9,pd12,pd23)
    
    data = cbind(data9,data12,data23)
    rownames(data) = rownames(data11)
    
    ##创建
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK11D_c9c12c23_monocle.rds')
    write.table(deg, file = "CK11D_c9c12c23_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  
  if(TRUE){
    fd <- as.data.frame(rownames(data11))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd4,pd7,pd8,pd13,pd15,pd21,pd24)
    
    data = cbind(data4,data7,data8,data13,data15,data21,data24)
    rownames(data) = rownames(data11)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK11D_c4c7c8c13c15c21c24_monocle.rds')
    write.table(deg, file = "CK11D_c4c7c8c13c15c21c24_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
  
  if(TRUE){
    fd <- as.data.frame(rownames(data11))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd17,pd18,pd20,pd25)
    
    data = cbind(data17,data18,data20,data25)
    rownames(data) = rownames(data11)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK11D_c17c18c20c25_monocle.rds')
    write.table(deg, file = "CK11D_c17c18c20c25_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}

# MCE
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data18 <- readRDS('CK3D_cluster18_count.rds')
  pd18 = metadata[colnames(data18),]
  pd18 = as.data.frame(pd18)
  pd18$celltype = 'CK3D_MCE'
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data16 <- readRDS('CK5D_cluster16_count.rds')
  pd16 = metadata[colnames(data16),]
  pd16 = as.data.frame(pd16)
  pd16$celltype = 'CK5D_MCE'
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data19 <- readRDS('CK7D_cluster19_count.rds')
  pd19 = metadata[colnames(data19),]
  pd19 = as.data.frame(pd19)
  pd19$celltype = 'CK7D_MCE'
  
  fd <- as.data.frame(rownames(data16))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd18,pd16,pd19)
  
  data = cbind(data18,data16,data19)
  rownames(data) = rownames(data16)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_MCE_3D5D7D_monocle.rds')
  write.table(deg, file = "CK_MCE_3D5D7D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  if(TRUE){
    fd <- as.data.frame(rownames(data16))
    colnames(fd) = 'GeneName'
    rownames(fd) = fd$GeneName
    
    pd = rbind(pd16,pd19)
    
    data = cbind(data16,data19)
    rownames(data) = rownames(data16)
    
    pd <- new("AnnotatedDataFrame", data = pd)
    fd <- new("AnnotatedDataFrame", data = fd)
    
    cds <- newCellDataSet(data,
                          phenoData = pd, 
                          featureData = fd,
                          expressionFamily = negbinomial.size())
    
    cds = estimateSizeFactors(cds)
    cds = estimateDispersions(cds)
    
    cds = detectGenes(cds, min_expr = 0.1)
    expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
    
    diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
    
    deg = subset(diff, qval < 0.01)
    deg = deg[order(deg$qval, decreasing = FALSE),]
    
    ordergene <- rownames(deg[1:2000,])
    cds = setOrderingFilter(cds, ordergene)
    
    cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
    cds = orderCells(cds)
    
    saveRDS(cds, 'CK_MCE_5D7D_monocle.rds')
    write.table(deg, file = "CK_MCE_5D7D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    
  }
    
}

# CZE
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data20 <- readRDS('CK3D_cluster20_count.rds')
  pd20 = metadata[colnames(data20),]
  pd20 = as.data.frame(pd20)
  pd20$celltype = 'CK3D_CZE'
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data14 <- readRDS('CK5D_cluster14_count.rds')
  pd14 = metadata[colnames(data14),]
  pd14 = as.data.frame(pd14)
  pd14$celltype = 'CK5D_CZE'
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data16 <- readRDS('CK7D_cluster16_count.rds')
  pd16 = metadata[colnames(data16),]
  pd16 = as.data.frame(pd16)
  pd16$celltype = 'CK7D_CZE'
  
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  data2 <- readRDS('CK9D_cluster2_count.rds')
  pd2 = metadata[colnames(data2),]
  pd2 = as.data.frame(pd2)
  pd2$celltype = 'CK9D_CZE'
  
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  data9 <- readRDS('CK11D_cluster9_count.rds')
  pd9 = metadata[colnames(data9),]
  pd9 = as.data.frame(pd9)
  pd9$celltype = 'CK11D_CZE'
  
  fd <- as.data.frame(rownames(data16))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd20,pd14,pd16,pd2,pd9)
  
  data = cbind(data20,data14,data16,data2,data9)
  rownames(data) = rownames(data16)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_CZE_3D5D7D9D11D_monocle.rds')
  write.table(deg, file = "CK_CZE_3D5D7D9D11D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

}

# PEN
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data18 <- readRDS('CK3D_cluster18_count.rds')
  pd18 = metadata[colnames(data18),]
  pd18 = as.data.frame(pd18)
  pd18$celltype = 'CK3D_PEN_c18'
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data3 <- readRDS('CK5D_cluster3_count.rds')
  pd3 = metadata[colnames(data3),]
  pd3 = as.data.frame(pd3)
  pd3$celltype = 'CK5D_PEN_c3'
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data4 <- readRDS('CK7D_cluster3_count.rds')
  pd4 = metadata[colnames(data4),]
  pd4 = as.data.frame(pd4)
  pd4$celltype = 'CK7D_PEN_c3'
  
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  data1 <- readRDS('CK9D_cluster1_count.rds')
  pd1 = getPD(data = data1, celltype = 'CK9D_c1')
  data9 <- readRDS('CK9D_cluster9_count.rds')
  pd9 = getPD(data = data9, celltype = 'CK9D_c9')
  data11 <- readRDS('CK9D_cluster11_count.rds')
  pd11 = getPD(data = data11, celltype = 'CK9D_c11')
  data15 <- readRDS('CK9D_cluster15_count.rds')
  pd15 = getPD(data = data15, celltype = 'CK9D_c15')
  data16 <- readRDS('CK9D_cluster16_count.rds')
  pd16 = getPD(data = data16, celltype = 'CK9D_c16')
  data19 <- readRDS('CK9D_cluster19_count.rds')
  pd19 = getPD(data = data19, celltype = 'CK9D_c19')
  data20 <- readRDS('CK9D_cluster20_count.rds')
  pd20 = getPD(data = data20, celltype = 'CK9D_c20')
  
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  data1ck11d <- readRDS('CK11D_cluster1_count.rds')
  data2ck11d <- readRDS('CK11D_cluster2_count.rds')
  data6ck11d <- readRDS('CK11D_cluster6_count.rds')
  data16ck11d <- readRDS('CK11D_cluster16_count.rds')
  data19ck11d <- readRDS('CK11D_cluster19_count.rds')
  data22ck11d <- readRDS('CK11D_cluster22_count.rds')
  pd1ck11d = getPD(data = data1ck11d, celltype = 'CK11D_c1')
  pd2ck11d = getPD(data = data2ck11d, celltype = 'CK11D_c2')
  pd6ck11d = getPD(data = data6ck11d, celltype = 'CK11D_c6')
  pd16ck11d = getPD(data = data16ck11d, celltype = 'CK11D_c16')
  pd19ck11d = getPD(data = data19ck11d, celltype = 'CK11D_c19')
  pd22ck11d = getPD(data = data22ck11d, celltype = 'CK11D_c22')

  fd <- as.data.frame(rownames(data16))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd18,pd3,pd4,pd1,pd9,pd11,pd15,pd16,pd19,pd20,pd1ck11d,pd2ck11d,pd6ck11d,pd16ck11d,pd19ck11d,pd22ck11d)
  
  data = cbind(data18,data3,data4,data1,data9,data11,data15,data16,data19,data20,data1ck11d,data2ck11d,data6ck11d,data16ck11d,data19ck11d,data22ck11d)
  rownames(data) = rownames(data16)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_PEN_3D5D7D9D11D_monocle.rds')
  write.table(deg, file = "CK_PEN_3D5D7D9D11D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

# CZSC-phloem-xylem
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data19 <- readRDS('CK3D_cluster19_count.rds')
  pd19 = getPD(data = data19, celltype = 'CK3D_CZSC_c19')
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data17 <- readRDS('CK5D_cluster17_count.rds')
  pd17 = getPD(data = data17, celltype = 'CK5D_CZSC_c17')
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data22 <- readRDS('CK7D_cluster22_count.rds')
  pd22 = getPD(data = data22, celltype = 'CK7D_CZSC_c22')
  
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  data20 <- readRDS('CK9D_cluster22_count.rds')
  pd20 = getPD(data = data20, celltype = 'CK9D_CZSC_c22')
  
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  data18 <- readRDS('CK11D_cluster18_count.rds')
  pd18 = getPD(data = data18, celltype = 'CK11D_CZSC_c18')
  
  fd <- as.data.frame(rownames(data18))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd19,pd17,pd22,pd20,pd18)
  
  data = cbind(data19,data17,data22,data20,data18)
  rownames(data) = rownames(data18)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_CZSC_3D5D7D9D11D_monocle.rds')
  write.table(deg, file = "CK_CZSC_3D5D7D9D11D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

# SC-oi
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data19 <- readRDS('CK3D_cluster1_count.rds')
  pd19 = getPD(data = data19, celltype = 'CK3D_CZSC_c1')
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data17 <- readRDS('CK5D_cluster1_count.rds')
  pd17 = getPD(data = data17, celltype = 'CK5D_CZSC_c1')
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data22 <- readRDS('CK7D_cluster1_count.rds')
  pd22 = getPD(data = data22, celltype = 'CK7D_CZSC_c1')
  
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  data20 <- readRDS('CK9D_cluster13_count.rds')
  pd20 = getPD(data = data20, celltype = 'CK9D_CZSC_c13')
  
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  data18 <- readRDS('CK11D_cluster14_count.rds')
  pd18 = getPD(data = data18, celltype = 'CK11D_CZSC_c14')
  
  fd <- as.data.frame(rownames(data18))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd19,pd17,pd22,pd20,pd18)
  
  data = cbind(data19,data17,data22,data20,data18)
  rownames(data) = rownames(data18)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_SCoi_3D5D7D9D11D_monocle.rds')
  write.table(deg, file = "CK_SCoi_3D5D7D9D11D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

# EP
if(TRUE){
  metadata = read.table("CK3D_metadata_res1.txt", header = TRUE)
  data19 <- readRDS('CK3D_cluster3_count.rds')
  pd19 = getPD(data = data19, celltype = 'CK3D_c3')
  
  metadata = read.table("CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data17 <- readRDS('CK5D_cluster11_count.rds')
  pd17 = getPD(data = data17, celltype = 'CK5D_c11')
  
  metadata = read.table("CK7D_metadata_res1.txt", header = TRUE)
  data11 <- readRDS('CK7D_cluster11_count.rds')
  pd11 = getPD(data = data11, celltype = 'CK7D_c11')
  data17 <- readRDS('CK7D_cluster17_count.rds')
  pd17 = getPD(data = data17, celltype = 'CK7D_c17')
  data25 <- readRDS('CK7D_cluster25_count.rds')
  pd25 = getPD(data = data25, celltype = 'CK7D_c25')
  
  metadata = read.table("CK9D_metadata_res1.txt", header = TRUE)
  data10 <- readRDS('CK9D_cluster10_count.rds')
  pd10 = getPD(data = data10, celltype = 'CK9D_c10')
  data21 <- readRDS('CK9D_cluster21_count.rds')
  pd21 = getPD(data = data21, celltype = 'CK9D_c21')
  
  metadata = read.table("CK11D_metadata_res1_subcluster.txt", header = TRUE)
  data4CK11D <- readRDS('CK11D_cluster4_count.rds')
  data7CK11D <- readRDS('CK11D_cluster7_count.rds')
  data8CK11D <- readRDS('CK11D_cluster8_count.rds')
  data13CK11D <- readRDS('CK11D_cluster13_count.rds')
  data15CK11D <- readRDS('CK11D_cluster15_count.rds')
  data21CK11D <- readRDS('CK11D_cluster21_count.rds')
  data24CK11D <- readRDS('CK11D_cluster24_count.rds')
  pd4CK11D = getPD(data = data4CK11D, celltype = 'CK11D_c4')
  pd7CK11D = getPD(data = data7CK11D, celltype = 'CK11D_c7')
  pd8CK11D = getPD(data = data8CK11D, celltype = 'CK11D_c8')
  pd13CK11D = getPD(data = data13CK11D, celltype = 'CK11D_c13')
  pd15CK11D = getPD(data = data15CK11D, celltype = 'CK11D_c15')
  pd21CK11D = getPD(data = data21CK11D, celltype = 'CK11D_c21')
  pd24CK11D = getPD(data = data24CK11D, celltype = 'CK11D_c24')

  fd <- as.data.frame(rownames(data11))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd19,pd17,pd11,pd17,pd25,pd10,pd21,pd4CK11D,pd7CK11D,pd8CK11D,pd13CK11D,pd15CK11D,pd21CK11D,pd24CK11D)
  
  data = cbind(data19,data17,data11,data17,data25,data10,data21,data4CK11D,data7CK11D,data8CK11D,data13CK11D,data15CK11D,data21CK11D,data24CK11D)
  rownames(data) = rownames(data11)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_EP_3D5D7D9D11D_monocle.rds')
  write.table(deg, file = "CK_EP_3D5D7D9D11D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

# SUS
if(TRUE){
  metadata = read.table("CK3D_bchi/CK3D_metadata_res1.txt", header = TRUE)
  data1 <- readRDS('CK3D_bchi/CK3D_cluster16_count.rds')
  pd1 = getPD(data = data1, celltype = 'CK3D_c16')
  
  metadata = read.table("CK5D_bchi/CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data2 <- readRDS('CK5D_bchi/CK5D_cluster15_count.rds')
  pd2 = getPD(data = data2, celltype = 'CK5D_c15')
  
  metadata = read.table("CK7D_bchi/CK7D_metadata_res1.txt", header = TRUE)
  data3 <- readRDS('CK7D_bchi/CK7D_cluster21_count.rds')
  pd3 = getPD(data = data3, celltype = 'CK7D_c21')
  
  fd <- as.data.frame(rownames(data1))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd1,pd2,pd3)
  
  data = cbind(data1,data2,data3)
  rownames(data) = rownames(data1)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_SUS_3D5D7D_monocle.rds')
  write.table(deg, file = "CK_SUS_3D5D7D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}

# SUS + CK9D SC
if(TRUE){
  metadata = read.table("CK3D_bchi/CK3D_metadata_res1.txt", header = TRUE)
  data1 <- readRDS('CK3D_bchi/CK3D_cluster16_count.rds')
  pd1 = getPD(data = data1, celltype = 'CK3D_SUS_c16')
  
  metadata = read.table("CK5D_bchi/CK5D_metadata_res1_subcluster.txt", header = TRUE)
  data2 <- readRDS('CK5D_bchi/CK5D_cluster15_count.rds')
  pd2 = getPD(data = data2, celltype = 'CK5D_SUS_c15')
  
  metadata = read.table("CK7D_bchi/CK7D_metadata_res1.txt", header = TRUE)
  data3 <- readRDS('CK7D_bchi/CK7D_cluster21_count.rds')
  pd3 = getPD(data = data3, celltype = 'CK7D_SUS_c21')
  
  metadata = read.table("CK9D_bchi/CK9D_metadata_res1.txt", header = TRUE)
  data4 <- readRDS('CK9D_bchi/CK9D_cluster7_count.rds')
  pd4 = getPD(data = data4, celltype = 'CK9D_SC_c7')
  
  fd <- as.data.frame(rownames(data1))
  colnames(fd) = 'GeneName'
  rownames(fd) = fd$GeneName
  
  pd = rbind(pd1,pd2,pd3,pd4)
  
  data = cbind(data1,data2,data3,data4)
  rownames(data) = rownames(data1)
  
  pd <- new("AnnotatedDataFrame", data = pd)
  fd <- new("AnnotatedDataFrame", data = fd)
  
  cds <- newCellDataSet(data,
                        phenoData = pd, 
                        featureData = fd,
                        expressionFamily = negbinomial.size())
  
  cds = estimateSizeFactors(cds)
  cds = estimateDispersions(cds)
  
  cds = detectGenes(cds, min_expr = 0.1)
  expressed_genes = rownames(subset(fData(cds), num_cells_expressed >=10))
  
  diff = differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~celltype", cores = 1)
  
  deg = subset(diff, qval < 0.01)
  deg = deg[order(deg$qval, decreasing = FALSE),]
  
  ordergene <- rownames(deg[1:2000,])
  cds = setOrderingFilter(cds, ordergene)
  
  cds = reduceDimension(cds, max_components = 2, method = 'DDRTree')
  cds = orderCells(cds)
  
  saveRDS(cds, 'CK_SUS_3D5D7D9D_monocle.rds')
  write.table(deg, file = "CK_SUS_3D5D7D9D_monocle_DEG.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  
}
