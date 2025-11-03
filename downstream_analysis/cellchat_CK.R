library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)

setwd('scripts/')

# CellChatDB
if(TRUE){
  LR_pair_bchi = read.table('LR_pair_atha_bchi.txt', sep = '\t', header = TRUE)
  LR_pair_bchi = LR_pair_bchi[LR_pair_bchi$source != 'orthologs',]
  LR_pair_bchi = LR_pair_bchi[LR_pair_bchi$Ligands_bchi != '',]
  LR_pair_bchi = LR_pair_bchi[LR_pair_bchi$Receptors_bchi != '',]
  LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == 'None','Ligands_atha_symbol'] = LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == 'None','Ligands']
  LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == '','Ligands_atha_symbol'] = LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == '','Ligands']
  LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == 'None','Receptors_atha_symbol'] = LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == 'None','Receptors']
  LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == '','Receptors_atha_symbol'] = LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == '','Receptors']
  LR_pair_bchi <- LR_pair_bchi %>% distinct(Ligands, Receptors, .keep_all = TRUE)
  head(LR_pair_bchi)
  table(LR_pair_bchi$source)
  
  db = list()
  db$interaction = data.frame(
    interaction_name = paste(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi, sep = '->'),
    pathway_name = paste(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol, sep = '->'),
    ligand = LR_pair_bchi$Ligands_bchi,
    receptor = LR_pair_bchi$Receptors_bchi,
    agonist = '',
    antagonist = '',
    co_A_receptor = '',
    co_I_receptor = '',
    evidence = LR_pair_bchi$source,
    annotation = '',
    interaction_name_2 = paste(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol, sep = '->')
  )
  rownames(db$interaction) = paste(LR_pair_bchi$Ligands, LR_pair_bchi$Receptors, sep = '->')
  
  db$cofactor = data.frame(
    cofactor1 = '',
    cofactor2 = '',
    cofactor3 = '',
    cofactor4 = ''
  )
  
  db$complex = data.frame(
    subunit_1 = '',
    subunit_2 = '',
    subunit_3 = '',
    subunit_4 = ''
  )
  
  db$geneInfo = data.frame(
    Symbol = c(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi),
    Name = c(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol),
    EntrezGene.ID = '',
    Ensembl.Gene.ID = c(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi),
    MGI.ID = '',
    Gene.group.name = ''
  )
  showDatabaseCategory(db)
  
  unique(db$interaction$annotation)
  
}

# CK3D
if(TRUE){
  sce = readRDS("CK3D_merge_res1.rds")
  cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")
  
  CellChatDB.use <- db # CellChatDB.use <- subsetDB(db, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "../revised_data/CK3D_net_lr.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, "../revised_data/CK3D_net_pathway.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  saveRDS(object = cellchat, file = "../revised_data/cellchat_CK3D.rds")
}

# CK5D
if(TRUE){
  sce = readRDS("CK5D_merge_res1_subcluster.rds")
  cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")
  
  CellChatDB.use <- db # CellChatDB.use <- subsetDB(db, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "../revised_data/CK5D_net_lr.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, "../revised_data/CK5D_net_pathway.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  saveRDS(object = cellchat, file = "../revised_data/cellchat_CK5D.rds")
}

# CK7D
if(TRUE){
  sce = readRDS("CK7D_merge_res1.rds")
  cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")
  
  CellChatDB.use <- db # CellChatDB.use <- subsetDB(db, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "../revised_data/CK7D_net_lr.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, "../revised_data/CK7D_net_pathway.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  saveRDS(object = cellchat, file = "../revised_data/cellchat_CK7D.rds")
}

# CK9D
if(TRUE){
  sce = readRDS("CK9D_merge_res1.rds")
  cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")
  
  CellChatDB.use <- db # CellChatDB.use <- subsetDB(db, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat) 
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "../revised_data/CK9D_net_lr.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, "../revised_data/CK9D_net_pathway.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  saveRDS(object = cellchat, file = "../revised_data/cellchat_CK9D.rds")
}

# CK11D
if(TRUE){
  sce = readRDS("CK11D_merge_res1_subcluster.rds")
  cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")
  
  CellChatDB.use <- db # CellChatDB.use <- subsetDB(db, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use # set the used database in the object
  
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  df.net <- subsetCommunication(cellchat)
  write.csv(df.net, "../revised_data/CK11D_net_lr.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  cellchat <- computeCommunProbPathway(cellchat)
  df.netp <- subsetCommunication(cellchat, slot.name = "netP")
  write.csv(df.netp, "../revised_data/CK11D_net_pathway.txt", sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  saveRDS(object = cellchat, file = "../revised_data/cellchat_CK11D.rds")
}
