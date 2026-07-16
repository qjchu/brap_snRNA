###############################################################################
# Purpose: Perform pseudotime trajectory analysis and downstream visualisation
#          using a pre-computed Monocle CellDataSet (CDS). This example uses
#          endosperm clusters from TT time series (5D, 7D, 8D, 9D, 11D).
# Steps:
#   1. Load CDS object and generate trajectory plots (pseudotime, cell type, state).
#   2. Relabel cell types (PCD vs active) and replot.
#   3. Plot time-course density and facet plots.
#   4. Visualise key genes along pseudotime.
#   5. Perform differential expression test along pseudotime and generate heatmap.
#   6. Cluster heatmap genes and run GO enrichment for each cluster.
###############################################################################

library(monocle)
library(dplyr)
library(Seurat)
library(ggplot2)
library(svglite)
library(ggsci)
library(clusterProfiler)
library(org.Brapaoleracea.eg.db)

# Set working directory (adjust as needed)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# 1. Load pre-computed Monocle CDS object
# ----------------------------------------------------------------------------
cds <- readRDS("../data/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle.rds")

# ----------------------------------------------------------------------------
# 2. Trajectory plots: Pseudotime, State, and cell types
# ----------------------------------------------------------------------------
# Plot by Pseudotime
p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1,
                           show_backbone = TRUE, alpha = 0.5) +
  scale_color_bs5("red") +
  theme(legend.position = "right")
ggsave("../plots/monocle_time.png", p1, height = 4, width = 7, bg = "transparent")
ggsave("../plots/monocle_time.svg", p1, height = 4, width = 7, bg = "transparent")

# Plot by State
plot_cell_trajectory(cds, color_by = "State", size = 1,
                     show_backbone = TRUE, alpha = 0.5) +
  theme(legend.position = "right")

# Plot by original celltype (before relabelling)
plot_cell_trajectory(cds, color_by = "celltype", size = 1,
                     show_backbone = TRUE, alpha = 0.5) +
  theme(legend.position = "right")

# ----------------------------------------------------------------------------
# 3. Relabel cell types to simplify (e.g., PCD vs active endosperm)
# ----------------------------------------------------------------------------
# Create a new column 'celltype2' with more descriptive labels
cds@phenoData@data[["celltype2"]] <- cds@phenoData@data[["celltype"]]
cds@phenoData@data[["celltype2"]] <- gsub("TT11D_Endo_c17", "TT11D_Endosperm", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT5D_Endo_c18", "TT5D_Endosperm_PCD", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT5D_Endo_c24", "TT5D_Endosperm_active", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT7D_Endo_c11", "TT7D_Endosperm_PCD", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT7D_Endo_c25", "TT7D_Endosperm_active", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT8D_Endo_c15", "TT8D_Endosperm_PCD", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT8D_Endo_c19", "TT8D_Endosperm_PCD", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT8D_Endo_c22", "TT8D_Endosperm_active", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT8D_Endo_c24", "TT8D_Endosperm_active", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT9D_Endo_c20", "TT9D_Endosperm_PCD", .)
cds@phenoData@data[["celltype2"]] <- gsub("TT9D_Endo_c28", "TT9D_Endosperm_active", .)

# Further simplify to "Endosperm_PCD" and "Endosperm_active" (ignore time)
cds@phenoData@data[["celltype3"]] <- cds@phenoData@data[["celltype2"]]
cds@phenoData@data[["celltype3"]] <- gsub("TT\\d+D_Endosperm_.*", "Endosperm_PCD",
                                           gsub("TT\\d+D_Endosperm_active", "Endosperm_active", .))

# Plot by simplified celltype3
p_type <- plot_cell_trajectory(cds, color_by = "celltype3", size = 1,
                               show_backbone = TRUE, alpha = 0.5) +
  theme(legend.position = "right") +
  scale_color_igv()
ggsave("../plots/monocle_type.svg", p_type, height = 4, width = 7, bg = "transparent")

# ----------------------------------------------------------------------------
# 4. Time-course visualisation
# ----------------------------------------------------------------------------
# Extract time from celltype (e.g., "TT5D_..." -> "TT5D")
cds@phenoData@data$time <- gsub("_.*", "", cds@phenoData@data$celltype)
cds@phenoData@data$time <- factor(cds@phenoData@data$time,
                                  levels = c("TT5D", "TT7D", "TT8D", "TT9D", "TT11D"))

# Colour palette for time points
time_colours <- c("#BDB76B", "#FFD700", "#FFA500", "#DAA520", "#B8860B")

# Trajectory coloured by time
p_time <- plot_cell_trajectory(cds, color_by = "time", size = 1,
                               show_backbone = TRUE) +
  scale_color_manual(values = time_colours) +
  theme(legend.position = "right")
ggsave("../plots/monocle_time_colour.svg", p_time, height = 4, width = 7, bg = "transparent")

# Facet by time
p_facet <- plot_cell_trajectory(cds, color_by = "time", size = 1,
                                show_backbone = TRUE) +
  scale_color_manual(values = time_colours) +
  theme(legend.position = "right") +
  facet_wrap(~time, nrow = 1)
ggsave("../plots/monocle_facet.svg", p_facet, height = 4, width = 15, bg = "transparent")
ggsave("../plots/monocle_facet.png", p_facet, height = 4, width = 15, bg = "transparent")

# Density plot of Pseudotime per time point
p_density <- ggplot(pData(cds), aes(Pseudotime, colour = time, fill = time)) +
  geom_density(bw = 0.5, size = 1, alpha = 0.5) +
  scale_color_manual(values = time_colours) +
  scale_fill_manual(values = time_colours) +
  theme_classic2()
ggsave("../plots/monocle_density.svg", p_density, height = 4, width = 7, bg = "transparent")

# ----------------------------------------------------------------------------
# 5. Visualise key genes along pseudotime
# ----------------------------------------------------------------------------
# Example: State 3 markers (custom gene list)
keygenes <- c("Bch06G010520", "Bch07G028630", "Bch09G069080",
              "Bch01G050010", "Bch05G049660", "Bch02G027190", "Bch07G045250")
plot_genes_in_pseudotime(cds[keygenes, ], color_by = "celltype")
plot_genes_in_pseudotime(cds[keygenes, ], color_by = "State")
plot_genes_in_pseudotime(cds[keygenes, ], color_by = "time")
plot_genes_in_pseudotime(cds[keygenes, ], color_by = "Pseudotime")

# ----------------------------------------------------------------------------
# 6. Differential expression along pseudotime and heatmap
# ----------------------------------------------------------------------------
# Load pre-computed DEG table (or compute on the fly)
deg <- read.table("../data/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_DEG.txt", header = TRUE)

# Test genes that vary with pseudotime (using top 2000 DEGs)
time_diff <- differentialGeneTest(cds[rownames(deg[1:2000, ]), ],
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
time_diff <- time_diff[time_diff$qval < 0.05, c(5, 2, 3, 4, 1, 6, 7)]

# Generate heatmap with 5 gene clusters
p_heatmap <- plot_pseudotime_heatmap(cds[time_diff$GeneName, ],
                                     num_clusters = 5,
                                     cores = 1,
                                     show_rownames = FALSE,
                                     return_heatmap = TRUE)
ggsave("../plots/monocle_heatmap.svg", p_heatmap, height = 10, width = 4, bg = "transparent")
ggsave("../plots/monocle_heatmap.png", p_heatmap, height = 10, width = 4, bg = "transparent")

# ----------------------------------------------------------------------------
# 7. GO enrichment for each gene cluster (from heatmap)
# ----------------------------------------------------------------------------
# Extract gene clusters from the heatmap tree
clusters <- cutree(p_heatmap$tree_row, k = 5)
clustering <- data.frame(Gene_Clusters = clusters, gene = names(clusters))

# Loop over each cluster and run GO enrichment
for (i in 1:5) {
  genes <- rownames(clustering[clustering$Gene_Clusters == i, ])
  # Handle genes with "Bo" prefix (add ".1" if needed)
  genes[grepl("^Bo", genes)] <- paste0(genes[grepl("^Bo", genes)], ".1")
  
  ego <- enrichGO(gene = genes,
                  OrgDb = org.Brapaoleracea.eg.db,
                  ont = "ALL",
                  keyType = "GID",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.05)
  if (!is.null(ego) && nrow(ego@result) > 0) {
    write.table(as.data.frame(ego@result),
                file = paste0("../plots/monocle_GO_cluster", i, ".txt"),
                sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}

# Optional: visualise GO results (if needed)
# dotplot(ego, split = "ONTOLOGY") + facet_grid(ONTOLOGY~., scales = "free")
