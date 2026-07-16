###############################################################################
# Purpose: Build a custom CellChat ligand-receptor database from orthology
#          mapping, then run CellChat analysis on a single sample (CK3D) to
#          infer cell-cell communication.
# Steps:
#   1. Read the ligand-receptor pair file (Arabidopsis -> B. rapa orthologs).
#   2. Format as a CellChat-compatible database.
#   3. Load Seurat object, create CellChat object, and run the pipeline:
#         - identify overexpressed genes and interactions
#         - compute communication probabilities
#         - filter, export results, and save the CellChat object.
# Note: The original script repeated the same code for multiple samples.
#       Here we keep only one representative (CK3D) for clarity.
###############################################################################

library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)

# Set working directory (adjust to your path)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# 1. Build custom CellChat database from orthology mapping
# ----------------------------------------------------------------------------
# Read the ligand-receptor pair file (contains Arabidopsis IDs, B. rapa orthologs,
# and other metadata). Path is a placeholder – adjust accordingly.
LR_pair_bchi <- read.table("../path/to/data/LR_pair_atha_bchi.txt", sep = "\t", header = TRUE)

# Clean up the data
LR_pair_bchi <- LR_pair_bchi[LR_pair_bchi$source != "orthologs", ]          # keep only curated pairs
LR_pair_bchi <- LR_pair_bchi[LR_pair_bchi$Ligands_bchi != "", ]             # remove missing ligand
LR_pair_bchi <- LR_pair_bchi[LR_pair_bchi$Receptors_bchi != "", ]           # remove missing receptor

# Replace placeholder symbols with actual gene names
LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == "None", "Ligands_atha_symbol"] <- 
  LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == "None", "Ligands"]
LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == "", "Ligands_atha_symbol"] <- 
  LR_pair_bchi[LR_pair_bchi$Ligands_atha_symbol == "", "Ligands"]
LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == "None", "Receptors_atha_symbol"] <- 
  LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == "None", "Receptors"]
LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == "", "Receptors_atha_symbol"] <- 
  LR_pair_bchi[LR_pair_bchi$Receptors_atha_symbol == "", "Receptors"]

# Keep only unique ligand-receptor pairs (based on original Arabidopsis IDs)
LR_pair_bchi <- LR_pair_bchi %>% distinct(Ligands, Receptors, .keep_all = TRUE)

# Create the database list structure expected by CellChat
db <- list()

# interaction table
db$interaction <- data.frame(
  interaction_name = paste(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi, sep = "->"),
  pathway_name = paste(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol, sep = "->"),
  ligand = LR_pair_bchi$Ligands_bchi,
  receptor = LR_pair_bchi$Receptors_bchi,
  agonist = "",
  antagonist = "",
  co_A_receptor = "",
  co_I_receptor = "",
  evidence = LR_pair_bchi$source,
  annotation = "",
  interaction_name_2 = paste(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol, sep = "->")
)
rownames(db$interaction) <- paste(LR_pair_bchi$Ligands, LR_pair_bchi$Receptors, sep = "->")

# cofactor, complex, and geneInfo tables (mostly empty placeholders)
db$cofactor <- data.frame(cofactor1 = "", cofactor2 = "", cofactor3 = "", cofactor4 = "")
db$complex <- data.frame(subunit_1 = "", subunit_2 = "", subunit_3 = "", subunit_4 = "")
db$geneInfo <- data.frame(
  Symbol = c(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi),
  Name = c(LR_pair_bchi$Ligands_atha_symbol, LR_pair_bchi$Receptors_atha_symbol),
  EntrezGene.ID = "",
  Ensembl.Gene.ID = c(LR_pair_bchi$Ligands_bchi, LR_pair_bchi$Receptors_bchi),
  MGI.ID = "",
  Gene.group.name = ""
)

# Optionally visualise database categories (for debugging)
showDatabaseCategory(db)

# ----------------------------------------------------------------------------
# 2. Run CellChat on a single sample (CK3D) using the custom database
#    (The original script repeated this for CK5D, CK7D, CK9D, CK11D.
#     We keep only one example; the pattern is identical.)
# ----------------------------------------------------------------------------
# Load Seurat object (path placeholder)
sce <- readRDS("../path/to/data/CK3D_merge_res1.rds")

# Create CellChat object, grouping by harmony clusters
cellchat <- createCellChat(object = sce, group.by = "harmony_clusters_res1")

# Assign the custom database
cellchat@DB <- db

# Preprocessing: subset data and identify overexpressed genes/interactions
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract and export the ligand-receptor interaction table
df_net <- subsetCommunication(cellchat)
write.csv(df_net, "../path/to/output/CK3D_net_lr.txt", sep = "\t",
          quote = FALSE, row.names = TRUE, col.names = TRUE)

# Compute pathway-level communication and export
cellchat <- computeCommunProbPathway(cellchat)
df_netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df_netp, "../path/to/output/CK3D_net_pathway.txt", sep = "\t",
          quote = FALSE, row.names = TRUE, col.names = TRUE)

# Save the full CellChat object for later use
saveRDS(object = cellchat, file = "../path/to/output/cellchat_CK3D.rds")
