###############################################################################
# Purpose: Compare expression of hormone-related gene sets (ABA, IAA, GA, XTH)
#          between CK (wild-type) and TT (mutant) across developmental timepoints.
#          For each cell type of interest (e.g., SC-SUS, PEN, CZSC), the script:
#            - Extracts average expression data from selected clusters.
#            - Maps Arabidopsis gene IDs to B. rapa orthologs.
#            - Generates boxplot + line plots showing expression trends over time.
#          Plots are saved as SVG files.
###############################################################################

library(RColorBrewer)
library(reshape2)
library(Mfuzz)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)

sessionInfo()

# Set working directory (adjust to your path)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# 1. Load average expression matrices for all timepoints (CK and TT)
#    These files contain gene expression averaged per cluster.
#    Paths are placeholders – adjust to your data location.
# ----------------------------------------------------------------------------
data_dir <- "../path/to/data/"

# CK timepoints
rna_averages0 <- read.table(paste0(data_dir, "CK3D/CK3D_RNA_harmony_clusters_res1_avg.txt"))
colnames(rna_averages0) <- paste0("CK3D_", colnames(rna_averages0))
rna_averages1 <- read.table(paste0(data_dir, "CK5D_subcluster/CK5D_subcluster_RNA_harmony_clusters_res1_avg.txt"))
colnames(rna_averages1) <- paste0("CK5D_", colnames(rna_averages1))
rna_averages2 <- read.table(paste0(data_dir, "CK7D/CK7D_RNA_harmony_clusters_res1_avg.txt"))
colnames(rna_averages2) <- paste0("CK7D_", colnames(rna_averages2))
rna_averages3 <- read.table(paste0(data_dir, "CK9D/CK9D_RNA_harmony_clusters_res1_avg.txt"))
colnames(rna_averages3) <- paste0("CK9D_", colnames(rna_averages3))
rna_averages4 <- read.table(paste0(data_dir, "CK11D_subcluster/CK11D_subcluster_RNA_harmony_clusters_res1_avg.txt"))
colnames(rna_averages4) <- paste0("CK11D_", colnames(rna_averages4))

# TT timepoints
rna_averages5 <- read.table(paste0(data_dir, "TT5D/TT5D_combined_RNA_clusters_res1_Louvain_avg.txt"))
colnames(rna_averages5) <- paste0("TT5D_", colnames(rna_averages5))
rna_averages6 <- read.table(paste0(data_dir, "TT7D/TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt"))
colnames(rna_averages6) <- paste0("TT7D_", colnames(rna_averages6))
rna_averages7 <- read.table(paste0(data_dir, "TT8D/TT8D_combined_RNA_harmony_clusters_res1_avg_rep1_rep3.txt"))
colnames(rna_averages7) <- paste0("TT8D_", colnames(rna_averages7))
rna_averages8 <- read.table(paste0(data_dir, "TT9D/TT9D_combined_RNA_harmony_clusters_res1_Louvain_avg.txt"))
colnames(rna_averages8) <- paste0("TT9D_", colnames(rna_averages8))
rna_averages9 <- read.table(paste0(data_dir, "TT11D/TT11D_combined_RNA_clusters_res1_avg.txt"))
colnames(rna_averages9) <- paste0("TT11D_", colnames(rna_averages9))

# ----------------------------------------------------------------------------
# 2. Helper function: Find B. rapa orthologs of Arabidopsis genes
#    Input:  atha_gene data.frame with Arabidopsis IDs and gene names/types.
#            exp_data: expression matrix (genes x samples).
#    Output: melted data frame with Expression, Time, Sample, and Gene labels.
# ----------------------------------------------------------------------------
FindBchiHomo <- function(atha_gene, exp_data) {
  # Load orthology mapping file (B. rapa -> A. thaliana)
  transfer <- read.table("../path/to/mapping/bchi_to_atha_besthit.txt", sep = " ")
  colnames(transfer) <- c("bchi_id", "atha_id")
  rownames(transfer) <- transfer$bchi_id
  
  # Keep only orthologs present in the Arabidopsis gene list
  bchi_gene <- transfer[transfer$atha_id %in% atha_gene$atha_gene, ]
  bchi_gene$symbol <- atha_gene[bchi_gene$atha_id, 2]   # add symbol/type column
  
  # Further filter to genes present in the expression matrix
  bchi_gene <- bchi_gene[bchi_gene$bchi_id %in% rownames(exp_data), ]
  data_subset <- exp_data[bchi_gene$bchi_id, ]
  data_subset <- as.data.frame(data_subset)
  data_subset$gene <- rownames(data_subset)
  
  # Melt to long format
  pdata <- reshape2::melt(data_subset)
  colnames(pdata) <- c("Gene", "Time", "Exp")
  pdata$sample <- gsub("_.*", "", pdata$Time)   # extract timepoint label (e.g., CK3D)
  
  # Create a label column (bchi_id + symbol) for potential annotation
  pdata$Label <- c(rep(NA, nrow(pdata) - length(bchi_gene$symbol)),
                   paste0(bchi_gene$bchi_id, "_", bchi_gene$symbol))
  return(pdata)
}

# ----------------------------------------------------------------------------
# 3. Define gene sets of interest (ABA, IAA, GA, XTH)
#    Each set is defined as a data.frame with Arabidopsis IDs and a type/name.
#    (Only one example is shown below; in the original script, multiple sets
#     are defined sequentially. They are all used later.)
# ----------------------------------------------------------------------------
# Example: ABA synthesis genes (from DeepSeek)
atha_gene_ABA_syn <- data.frame(
  atha_gene = c("AT5G67030", "AT3G14440", "AT1G30100", "AT3G24220", "AT1G78390",
                "AT1G52340", "AT1G16540", "AT2G27150", "AT1G60680", "AT4G34000",
                "AT3G19290"),
  type = c("ABA1", "NCED3", "NCED5", "NCED6", "NCED9", "ABA2", "ABA3", "AAO",
           "LOS5/ABA4", "ABF3", "ABF4")
)
rownames(atha_gene_ABA_syn) <- atha_gene_ABA_syn$atha_gene

# (Other gene sets – ABA degradation, IAA synthesis/degradation, GA synthesis/degradation,
#  XTH – are defined similarly in the original script; we keep them for completeness.)

# ----------------------------------------------------------------------------
# 4. For each cell type of interest, extract the relevant columns from the
#    average expression matrices, then generate plots for each gene set.
#    The original script repeats this for many cell types (SC-SUS, PEN, CZSC, etc.)
#    and saves plots to different output files.
#    For brevity, we show the structure for one cell type (CZSC) and note that
#    the same pattern applies to others.
# ----------------------------------------------------------------------------

# 4a. Example: CZSC cell type
#     CK CZSC clusters: CK3D c5, CK5D c6, CK7D c8, CK9D c4, CK11D c17
#     TT CZSC clusters: TT5D c8, TT7D c19, TT8D c2, TT9D c15, TT11D c1
if (TRUE) {
  # Build CK expression matrix for CZSC (columns selected by index)
  a <- cbind(rna_averages0[rownames(rna_averages0), 5],
             rna_averages1[rownames(rna_averages0), 6],
             rna_averages2[rownames(rna_averages0), 8],
             rna_averages3[rownames(rna_averages0), 4],
             rna_averages4[rownames(rna_averages0), 17])
  colnames(a) <- c("CK3D", "CK5D", "CK7D", "CK9D", "CK11D")
  rownames(a) <- rownames(rna_averages0)
  a <- na.omit(a)
  
  # Build TT expression matrix for CZSC
  b <- cbind(rna_averages5[rownames(rna_averages0), 8],
             rna_averages6[rownames(rna_averages0), 19],
             rna_averages7[rownames(rna_averages0), 2],
             rna_averages8[rownames(rna_averages0), 15],
             rna_averages8[rownames(rna_averages0), 1])   # Note: TT11D uses c1
  colnames(b) <- c("TT5D", "TT7D", "TT8D", "TT9D", "TT11D")
  rownames(b) <- rownames(rna_averages0)
  b <- na.omit(b)
  
  # For each gene set, generate plot
  # (Here we show ABA synthesis; other sets would be similar)
  pdata1 <- FindBchiHomo(atha_gene = atha_gene_ABA_syn, exp_data = a)
  pdata2 <- FindBchiHomo(atha_gene = atha_gene_ABA_syn, exp_data = b)
  
  # Optional: Wilcoxon test between CK and TT (not used for plotting)
  # wilcox.test(pdata1$Exp, pdata2$Exp, paired = FALSE)
  
  pdata <- rbind(pdata1, pdata2)
  pdata$sample <- substr(pdata$sample, 1, 2)   # keep only 'CK' or 'TT'
  
  # Create the plot: boxplot + points + lines for each gene
  p <- ggplot(pdata, aes(x = Time, y = Exp, color = Time)) +
    geom_boxplot(width = 0.5, cex = 1) +
    geom_point(alpha = 0.6) +
    geom_line(aes(group = Gene), linetype = 1, alpha = 0.6) +
    scale_color_manual(values = rep(c("#9ACD32", "#6B8E23", "#556B2F",
                                      "#32CD32", "#228B22", "#BDB76B",
                                      "#FFD700", "#FFA500", "#DAA520",
                                      "#B8860B"), 50)) +
    facet_grid(~ sample, space = "free", scales = "free") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = "darkgray"),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          legend.position = "none")
  
  # Save the plot (adjust path as needed)
  ggsave("../plots/CK_TT_CZSC_ABAsyn.svg", p, height = 3, width = 3.5, bg = "transparent")
  
  # (The original script saves plots for many other gene sets and cell types.
  #   The pattern is identical – we omit the repetitive code here for clarity.)
}

# ----------------------------------------------------------------------------
# 5. Additional cell types (SC-SUS, PEN, EP, SC-ii, SC-oi, CZSC-phloem-xylem, etc.)
#    are handled in the same way in the original script. The only difference
#    is the column indices chosen for each cell type and the output filenames.
#    The code structure is identical; we do not reproduce all blocks here.
#    The script also includes a version with geom_text_repel for labeling.
# ----------------------------------------------------------------------------

# End of script
