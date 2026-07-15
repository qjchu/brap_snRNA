###############################################################################
# Purpose: Generate violin plots (flipped) for nFeature_RNA and nCount_RNA
#          for each sample in datasets, with quantile summaries.
# Input:  Seurat object (RDS) containing merged samples.
# Output: SVG plots and summary statistics printed to console.
###############################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(svglite)

# ----------------------------------------------------------------------------
# 1. Load the Seurat object (adjust path to your data location)
# ----------------------------------------------------------------------------
# Path placeholder – replace with your actual relative or absolute path
rds_path <- "../data/CK_batch3_10samples/10samples/CK_10samples_batch3_metadata_res1.rds"
rds <- readRDS(rds_path)

# ----------------------------------------------------------------------------
# 2. Define sample grouping and factor levels (order for plotting)
# ----------------------------------------------------------------------------
# Extract sample name from orig.ident (e.g., "CK3D_1" -> "CK3D")
rds@meta.data$sample <- gsub(pattern = "_.*", replacement = "", x = rds@meta.data$orig.ident)
# Set desired order (reverse for flipped coordinates)
sample_levels <- rev(c("CK3D", "CK5D", "CK7D", "CK9D", "CK11D"))
rds@meta.data$sample <- factor(rds@meta.data$sample, levels = sample_levels)

# ----------------------------------------------------------------------------
# 3. Helper function: generate VlnPlot data, compute stats, and save plot
# ----------------------------------------------------------------------------
plot_QC_violin <- function(feature, ylim_max = NULL, colors, output_file) {
  
  # Generate the VlnPlot and extract the underlying data
  vp <- VlnPlot(rds, features = feature, group.by = "sample", pt.size = 0, combine = FALSE)
  plot_data <- vp[[1]]$data   # data frame with columns: "ident" (sample) and feature
  
  # Print quantiles and cell counts per sample
  cat("\n=== Summary for", feature, "===\n")
  for (s in sample_levels) {
    subset_data <- plot_data[plot_data$ident == s, feature]
    q <- quantile(subset_data, probs = c(0, 0.25, 0.5, 0.75, 1))
    cat(s, ":\n")
    print(q)
    cat("  n =", length(subset_data), "\n\n")
  }
  
  # Create the violin + boxplot (flipped)
  p <- ggplot(plot_data, aes(x = factor(ident), y = !!sym(feature), color = ident)) +
    geom_violin(scale = "width", adjust = 1, trim = TRUE, cex = 1.2) +
    geom_boxplot(width = 0.1, cex = 1.2) +
    scale_color_manual(values = colors) +
    coord_flip() +
    labs(x = NULL, y = NULL) +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      axis.text.x = element_text(colour = "black", size = 10, angle = 90, vjust = 0.5, hjust = 0.5),
      axis.line = element_line(colour = "grey60"),
      legend.position = "none"
    )
  
  if (!is.null(ylim_max)) {
    p <- p + ylim(0, ylim_max)
  }
  
  # Save the plot
  ggsave(output_file, plot = p, height = 4, width = 1.6)
  cat("Plot saved:", output_file, "\n")
}

# ----------------------------------------------------------------------------
# 4. Generate plots for nFeature_RNA and nCount_RNA
# ----------------------------------------------------------------------------
# Colour palette for CK samples (green shades)
ck_colors <- rev(c("#9ACD32", "#6B8E23", "#556B2F", "#32CD32", "#228B22"))

# nFeature_RNA
plot_QC_violin(
  feature = "nFeature_RNA",
  ylim_max = NULL,
  colors = ck_colors,
  output_file = "../data/CK_batch3_10samples/CK_vlnplot_nFeature_RNA.svg"
)

# nCount_RNA (with y‑axis limit)
plot_QC_violin(
  feature = "nCount_RNA",
  ylim_max = 10000,
  colors = ck_colors,
  output_file = "../data/CK_batch3_10samples/CK_vlnplot_nCount_RNA.svg"
)
