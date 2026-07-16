###############################################################################
# Purpose: Compare gene expression between TT (mutant) and CK (wild-type) for
#          matched cell types. Calculate log2 fold change (TT/CK) for each gene
#          in each cell type and generate volcano-style dot plots with gene labels.
# Steps:
#   1. Load average expression matrices for TT7D and CK7D.
#   2. For each cell type pair, compute logFC (genes with |logFC| > 2).
#   3. Combine results and visualise as a column/point plot (Figure 4A style).
###############################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(svglite)
library(ggpubr)
library(ggrepel)

# Set working directory (adjust to your path)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# 1. Helper function to compute DEGs between two matched cell types
# ----------------------------------------------------------------------------
# Input: expression matrices (genes x samples), column names for TT and CK,
#        and a cell type label.
# Output: data frame with logFC and metadata for genes with |logFC| > 2.
celltype_DEG <- function(TT_data, CK_data, TT_celltype, CK_celltype, celltype) {
  # Subset columns and combine
  temp_data <- cbind(TT_data[, TT_celltype, drop = FALSE],
                     CK_data[rownames(TT_data), CK_celltype, drop = FALSE])
  colnames(temp_data) <- c("TT", "CK")
  rownames(temp_data) <- rownames(TT_data)
  temp_data <- na.omit(as.data.frame(temp_data))
  
  # Keep genes with reasonable expression in both (sum > 0.5 and both > 0)
  temp_data <- temp_data[temp_data$TT + temp_data$CK > 0.5 &
                         temp_data$TT > 0 & temp_data$CK > 0, ]
  
  # Compute log2 fold change (TT / CK)
  temp_data$logFC <- log2(temp_data$TT / temp_data$CK)
  
  # Filter by |logFC| > 2
  temp_data <- temp_data[abs(temp_data$logFC) > 2, ]
  temp_data$cluster <- celltype
  temp_data$GeneID <- rownames(temp_data)
  
  return(temp_data)
}

# ----------------------------------------------------------------------------
# 2. Load expression matrices (average expression per cluster)
#    Paths are placeholders – adjust to your data location.
# ----------------------------------------------------------------------------
# CK7D (wild-type) – 25 clusters
CK7D <- read.table("../data/CK7D/CK7D_RNA_harmony_clusters_res1_avg.txt", header = TRUE)
colnames(CK7D) <- c("CK7D SC-oi (c1)", "CK7D SC (c2)", "CK7D PEN (c3)", "CK7D SC_Endosperm (c4)",
                    "CK7D SC-ii (c5)", "CK7D SC (c6)", "CK7D SC_Dividing_Cell (c7)",
                    "CK7D CZSC (c8)", "CK7D CZSC (c9)", "CK7D SC (c10)",
                    "CK7D EP_Dividing_Cell (c11)", "CK7D SC (c12)", "CK7D SC (c13)",
                    "CK7D SC-ii (c14)", "CK7D SC_Endosperm (c15)", "CK7D CZE (c16)",
                    "CK7D EP_Dividing_Cell (c17)", "CK7D SC-oi (c18)", "CK7D MCE (c19)",
                    "CK7D SC (c20)", "CK7D SC_SUS (c21)", "CK7D CZSC-phloem-xylem (c22)",
                    "CK7D SC (c23)", "CK7D CZSC (c24)", "CK7D EP (c25)")
# Remove genes with very low total expression
CK7D <- CK7D[rowSums(CK7D) > 1, ]

# TT7D (mutant) – 25 clusters
TT7D <- read.table("../data/TT7D/TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt", header = TRUE)
colnames(TT7D) <- c("TT7D SC-oi (c1)", "TT7D SC-ii (c2)", "TT7D SC (c3)",
                    "TT7D SC-endosperm (c4)", "TT7D SC (c5)", "TT7D SC-ii (c6)",
                    "TT7D SC (c7)", "TT7D SC (c8)", "TT7D SC-dividing-cell (c9)",
                    "TT7D SC-oi (c10)", "TT7D Endosperm-PCD (c11)", "TT7D SC (c12)",
                    "TT7D SC (c13)", "TT7D SC-ii (c14)", "TT7D SC (c15)",
                    "TT7D CZSC (c16)", "TT7D CZSC (c17)", "TT7D SC (c18)",
                    "TT7D CZSC (c19)", "TT7D SC-ii (c20)", "TT7D SC-oi (c21)",
                    "TT7D EP-dividing-cell (c22)", "TT7D CZSC-phloem-xylem (c23)",
                    "TT7D SC-SUS (c24)", "TT7D Endosperm-active (c25)")
TT7D <- TT7D[rowSums(TT7D) > 1, ]

# Ensure genes are in the same order
common_genes <- intersect(rownames(TT7D), rownames(CK7D))
TT7D <- TT7D[common_genes, ]
CK7D <- CK7D[common_genes, ]

# ----------------------------------------------------------------------------
# 3. Compute DEGs for each matched cell type pair
# ----------------------------------------------------------------------------
# Define TT and CK cell type names (must match column names exactly)
pairs <- list(
  list(TT = "TT7D Endosperm-active (c25)", CK = "CK7D PEN (c3)", label = "Endosperm-active"),
  list(TT = "TT7D CZSC (c16)", CK = "CK7D CZSC (c9)", label = "CZSC"),
  list(TT = "TT7D CZSC-phloem-xylem (c23)", CK = "CK7D CZSC-phloem-xylem (c22)", label = "CZSC-phloem-xylem"),
  list(TT = "TT7D EP-dividing-cell (c22)", CK = "CK7D EP_Dividing_Cell (c17)", label = "EP-dividing-cell"),
  list(TT = "TT7D Endosperm-PCD (c11)", CK = "CK7D PEN (c3)", label = "Endosperm-PCD"),
  list(TT = "TT7D SC-oi (c1)", CK = "CK7D SC-oi (c1)", label = "SC-oi"),
  list(TT = "TT7D SC-SUS (c24)", CK = "CK7D SC_SUS (c21)", label = "SC-SUS"),
  list(TT = "TT7D SC-dividing-cell (c9)", CK = "CK7D SC_Dividing_Cell (c7)", label = "SC_Dividing_Cell"),
  list(TT = "TT7D SC (c5)", CK = "CK7D SC (c2)", label = "SC1"),
  list(TT = "TT7D SC (c7)", CK = "CK7D SC (c10)", label = "SC2"),
  list(TT = "TT7D SC (c15)", CK = "CK7D SC (c12)", label = "SC3"),
  list(TT = "TT7D SC-ii (c2)", CK = "CK7D SC-ii (c5)", label = "SC-ii"),
  list(TT = "TT7D SC-endosperm (c4)", CK = "CK7D SC_Endosperm (c15)", label = "SC_Endosperm")
)

# Apply function and collect results
deg_list <- list()
for (p in pairs) {
  deg_list[[p$label]] <- celltype_DEG(TT_data = TT7D, CK_data = CK7D,
                                      TT_celltype = p$TT, CK_celltype = p$CK,
                                      celltype = p$label)
}
deg_all <- do.call(rbind, deg_list)

# ----------------------------------------------------------------------------
# 4. Visualisation (Figure 4A style)
# ----------------------------------------------------------------------------
# Prepare data for plotting: jitter x-coordinates and compute bar heights
deg_all$jittered_x <- jitter(as.numeric(factor(deg_all$cluster)), amount = 0.35)

# Define the order and labels for x-axis
celltypes <- c("Endosperm-active", "CZSC", "CZSC-phloem-xylem", "EP-dividing-cell",
               "Endosperm-PCD", "SC-oi", "SC-SUS", "SC_Dividing_Cell",
               "SC1", "SC2", "SC3", "SC-ii", "SC_Endosperm")

# Compute max/min logFC for each cell type (for background bars)
df_max <- data.frame(x = celltypes,
                     y = sapply(celltypes, function(ct) max(deg_all$logFC[deg_all$cluster == ct])))
df_min <- data.frame(x = celltypes,
                     y = sapply(celltypes, function(ct) min(deg_all$logFC[deg_all$cluster == ct])))
df_label <- data.frame(x = celltypes, y = 0, label = celltypes)

# Optionally annotate with gene symbols (load annotation file)
anno <- read.table("../data/besthit_anno.txt", sep = "\t", fill = TRUE, header = FALSE)
colnames(anno) <- c("bchi_gene", "atha_gene", "symbol", "des1", "des2")
rownames(anno) <- anno$bchi_gene
anno$symbol[anno$symbol == "None"] <- anno$bchi_gene[anno$symbol == "None"]

deg_all$symbol <- anno[deg_all$GeneID, "symbol"]
deg_all$symbol[is.na(deg_all$symbol)] <- deg_all$GeneID[is.na(deg_all$symbol)]

# Highlight specific genes of interest (custom list)
highlight_genes <- c("AGL62", "ICE1", "AGL91", "UBP13", "LEC1", "TT2",
                     "CYP78A5", "MYB65", "MET1", "Bch06G010590", "PAB8",
                     "XTH24", "HB21", "LCR84", "DRM2", "Bch04G013320",
                     "ENODL3", "LCR30")
deg_all$symbol2 <- NA
deg_all$symbol2[deg_all$symbol %in% highlight_genes] <- deg_all$symbol[deg_all$symbol %in% highlight_genes]

# Create the plot (final version, similar to fig4A_v2)
p <- ggplot() +
  # Background bars (max and min)
  geom_col(data = df_max, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  geom_col(data = df_min, aes(x = x, y = y), fill = "#dcdcdc", alpha = 0.6) +
  # Points for each gene
  geom_point(data = deg_all, aes(x = jittered_x, y = logFC, color = logFC), size = 2, alpha = 0.96) +
  # Gene labels (only for highlighted genes)
  geom_text_repel(data = deg_all, aes(x = jittered_x, y = logFC, label = symbol2),
                  force = 3, box.padding = 0.8,
                  arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
  # Tile labels at bottom
  geom_tile(data = df_label, aes(x = x, y = y), height = 3.6, color = "white",
            fill = "darkgreen", alpha = 0.6, show.legend = FALSE) +
  geom_text(data = df_label, aes(x = x, y = y, label = label),
            angle = 90, vjust = 0.5, hjust = 0.5, size = 5, color = "white") +
  scale_color_gradient2(low = "#00e600", high = "red") +
  labs(y = "Log2(FoldChange)") +
  theme_minimal() +
  ylim(-14, 18) +
  theme(axis.title = element_text(size = 13, color = "black"),
        axis.title.y = element_blank(),
        axis.line.x = element_line(color = "gray45", size = 1),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  coord_flip()

# Save the figure
ggsave("../plots/TT7D_VS_CK7D_celltype_fig4A.svg", p, height = 18, width = 8, bg = "transparent")
ggsave("../plots/TT7D_VS_CK7D_celltype_fig4A.png", p, height = 18, width = 8, bg = "transparent")

# ----------------------------------------------------------------------------
# 5. Export full DEG table for supplementary data
# ----------------------------------------------------------------------------
write.table(deg_all, file = "../tables/TT7D_VS_CK7D_celltype_DEG_supplement.txt",
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
