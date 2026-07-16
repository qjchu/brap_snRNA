###############################################################################
# Purpose: Create a network-style plot showing connections between cell clusters
#          (nodes) and transition relationships (curved links) for a single
#          experimental condition (here, CK time series).
#          The plot combines:
#            - Points representing clusters (size ~ cell proportion)
#            - Curved edges representing significant correlations/relationships
#              between clusters, with line width ~ correlation strength.
#          Data are pre‑computed in Excel sheets (cluster metadata and edge list).
# Input:  Two Excel files: "CK_clusters.xlsx" (node info) and
#         "CK_clusters_cor.xlsx" (edge info). Both must be in the same directory.
# Output: SVG file with the final plot.
###############################################################################

library(tidyverse)
library(readxl)
library(ggplot2)
library(ggpubr)
library(grid)
library(ggsci)
library(ggnewscale)   # for multiple colour scales

# Source external function for curved links (supplied in the same directory)
source("geom-curve2.R")

# ----------------------------------------------------------------------------
# 1. Read cluster metadata (nodes) for CK samples
#    (the script loads several metadata files but only the last one is used
#     to compute proportions; we keep a simplified example)
# ----------------------------------------------------------------------------
# Path placeholder – adjust to your actual data location
data_dir <- "../data/CK_batch3_10samples/"

# For demonstration, we read one metadata file (CK11D) and compute cluster proportions.
# In the original script, several files were loaded; the last one (CK11D_subcluster)
# was used for the proportions shown in the final plot.
metadata <- read.table(paste0(data_dir, "CK11D_subcluster/CK11D_metadata_res1_subcluster.txt"),
                       header = TRUE)
prop_table <- table(metadata$harmony_clusters_res1) / nrow(metadata)

# The original script then manually ordered clusters for each timepoint.
# We only keep the last ordering (CK11D) as an example.
# The manual ordering is used to set the node positions in the Excel file;
# we do not reproduce it here because the Excel files already contain the final
# coordinates (newX, newY) and cluster identities.

# ----------------------------------------------------------------------------
# 2. Read node and edge data from Excel files
# ----------------------------------------------------------------------------
# Node file: clusters with coordinates, cluster labels, and group (tissue type)
dat01 <- read_excel(paste0(data_dir, "CK_clusters.xlsx"))
dat01$newX <- dat01$newX - 0.5   # shift X coordinates for visual centering

# Edge file: relationships between clusters (startx, starty, endx, endy)
dat02 <- read_excel(paste0(data_dir, "CK_clusters_cor.xlsx"))
dat02$startx <- dat02$startx - 0.5
dat02$endx   <- dat02$endx   - 0.5

# Keep only the strongest edge per starting cluster (optional, as in original)
dat02 <- dat02 %>% arrange(desc(size2))
top_edges <- dat02 %>% group_by(startx, starty) %>% slice_head(n = 1) %>% ungroup()

# ----------------------------------------------------------------------------
# 3. Generate the plot
# ----------------------------------------------------------------------------
p <- ggplot() +
  # Draw curved edges (links)
  geom_curve2(data = top_edges, alpha = 0.5,
              aes(x = startx, xend = endx,
                  y = starty, yend = endy,
                  curvature = curvature,
                  size = size2,        # line width proportional to correlation strength
                  color = group),
              node.color = NA, ncp = 1, node.fill = NA) +
  scale_color_npg() +
  
  # New scale for node colours (points) – override the edge colour scale
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale("size") +
  
  # Draw cluster nodes as points (size ~ cell proportion/value3)
  geom_point(data = dat01,
             aes(x = newX, y = newY, size = value3, color = group)) +
  scale_size_continuous(range = c(0.5, 15)) +
  scale_linewidth_continuous(range = c(0, 0.1)) +
  
  # Set node colours by tissue group (Endosperm, EP, SC)
  scale_color_manual(values = c("Endosperm" = "#e64b35",
                                "EP"        = "#4dbbd5",
                                "SC"        = "#00a087")) +
  
  # Add cluster labels using annotate()
  annotate(geom = "label",
           x = dat01$newX, y = dat01$newY,
           label = dat01$celltype3,
           label.size = NA, fill = "grey95", size = 3, alpha = 0.5) +
  
  theme_void() +
  xlim(0, 5)

# Save the plot
ggsave("../plots/CK_clusters_link.svg", plot = p, height = 8, width = 12)
