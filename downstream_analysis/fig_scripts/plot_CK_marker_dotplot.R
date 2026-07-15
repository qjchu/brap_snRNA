###############################################################################
# Purpose: Generate a dotplot panel showing expression of conserved marker genes
#          across clusters from multiple CK (control) timepoints (CK3D, CK5D, CK7D, CK9D, CK11D).
# Workflow:
#   1. Load marker gene tables for each timepoint (from FindAllMarkers output).
#   2. For each cell type, intersect the top 100 marker genes across relevant clusters/timepoints
#      to identify conserved markers, then select top 3 per type.
#   3. Load Seurat objects for each timepoint and generate DotPlot data for the selected markers.
#   4. Combine data, add cell type annotations, and create a faceted dotplot (timepoint × cell type).
# Output: SVG file with the dotplot.
###############################################################################

library(dplyr)
library(Seurat)
library(ggplot2)
library(svglite)

# ----------------------------------------------------------------------------
# 1. Read marker gene tables for each CK timepoint
#    (assumes files are in a consistent directory structure)
# ----------------------------------------------------------------------------
# Path placeholders – adjust to your actual data location
data_dir <- "../data/CK_batch3_10samples/"

# Function to read and clean a marker file
read_marker_file <- function(file_path) {
  df <- read.table(file_path, header = TRUE, sep = "\t", quote = "")
  colnames(df) <- c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster',
                    'gene','atha_GeneID','atha_Symbol','computational_description','full_name')
  df <- df %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
  return(df)
}

# Read each timepoint
m1 <- read_marker_file(paste0(data_dir, "CK3D/CK3D_res1_harmony_cluster_markers_anno.txt"))
m2 <- read_marker_file(paste0(data_dir, "CK5D_subcluster/CK5D_res1_subcluster_harmony_cluster_markers_anno.txt"))
m3 <- read_marker_file(paste0(data_dir, "CK7D/CK7D_res1_harmony_cluster_markers_anno.txt"))
m4 <- read_marker_file(paste0(data_dir, "CK9D/CK9D_res1_harmony_cluster_markers_anno.txt"))
m5 <- read_marker_file(paste0(data_dir, "CK11D_subcluster/CK11D_res1_subcluster_harmony_cluster_markers_anno.txt"))

# Take top 100 markers per cluster for each timepoint
m1_top <- m1 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()
m2_top <- m2 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()
m3_top <- m3 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()
m4_top <- m4 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()
m5_top <- m5 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()

# ----------------------------------------------------------------------------
# 2. Define cell types and their corresponding cluster(s) across timepoints,
#    then intersect marker gene lists to get conserved markers.
#    For each type, we then manually select 3 representative genes.
# ----------------------------------------------------------------------------
# Define lists of cluster numbers for each cell type per timepoint
# (These values are specific to the CK dataset annotation)

# For brevity, we directly define the final marker gene sets (already curated)
markers <- c(
  # EP-dividing-cell
  "Bch01G042330", "Bch02G014800", "Bch06G006190",
  # EP
  "Bch01G045760", "Bch01G047930", "Bch07G013960",
  # CZE
  "Bch06G008080", "Bch06G009710", "Bch07G011060",
  # MCE
  "Bch03G067580", "Bch04G033350", "Bch05G008080",
  # PEN
  "Bch03G072300", "Bch04G014210", "Bch09G041620",
  # CZSC
  "Bch05G007030", "Bch01G000920", "Bch09G050570",
  # CZSC-phloem-xylem
  "Bch09G027690", "Bch02G031690", "Bch03G063320",
  # SC-ii
  "Bch01G013080", "Bch09G017230", "Bch07G022460",
  # SC-oi
  "Bch01G041200", "Bch07G015970", "Bch01G011550",
  # SC-SUS
  "Bch03G029400", "Bch05G001280", "Bch03G065950"
)

# Define cell type names in the same order (for annotation)
celltype_labels <- c(
  rep("EP-dividing-cell", 3),
  rep("EP", 3),
  rep("CZE", 3),
  rep("MCE", 3),
  rep("PEN", 3),
  rep("CZSC", 3),
  rep("CZSC-phloem-xylem", 3),
  rep("SC-ii", 3),
  rep("SC-oi", 3),
  rep("SC-SUS", 3)
)

# ----------------------------------------------------------------------------
# 3. Load Seurat objects for each CK timepoint and extract DotPlot data
# ----------------------------------------------------------------------------
rds_dir <- "../data/RDS/"   # adjust as needed

# List of timepoints and corresponding RDS files
timepoints <- c("CK3D", "CK5D", "CK7D", "CK9D", "CK11D")
rds_files <- c(
  "CK3D_merge_res1.rds",
  "CK5D_merge_res1_subcluster.rds",
  "CK7D_merge_res1.rds",
  "CK9D_merge_res1.rds",
  "CK11D_merge_res1_subcluster.rds"
)

# Function to extract DotPlot data from a Seurat object
get_dotplot_data <- function(rds_path, time_label) {
  sce <- readRDS(rds_path)
  dp <- DotPlot(sce, features = markers)
  data <- dp[["data"]]
  data$id <- paste0(time_label, "_", data$id)
  data$time <- time_label
  return(data)
}

# Collect all dotplot data into one data frame
dotplot_list <- list()
for (i in seq_along(timepoints)) {
  rds_path <- file.path(rds_dir, rds_files[i])
  dotplot_list[[i]] <- get_dotplot_data(rds_path, timepoints[i])
}
dotplot_data <- bind_rows(dotplot_list)

# ----------------------------------------------------------------------------
# 4. Tidy and annotate the combined data
# ----------------------------------------------------------------------------
# Order clusters manually (optional – keep the factor levels as in original script)
# Here we use a simplified version; the full list is long but can be adapted.

# For simplicity, we keep the default order (or we can set a custom order).
# The original script had a long vector of factor levels; we omit it here for brevity,
# but you can add it if needed. We'll just set levels to the unique ids in the data.
dotplot_data$id <- factor(dotplot_data$id, levels = sort(unique(dotplot_data$id)))

# Add cell type annotation based on marker genes
dotplot_data$celltype <- NA
for (i in seq_along(markers)) {
  dotplot_data$celltype[dotplot_data$features.plot == markers[i]] <- celltype_labels[i]
}
dotplot_data$celltype <- factor(dotplot_data$celltype,
                                levels = rev(c("EP-dividing-cell", "EP", "CZE", "MCE", "PEN",
                                               "CZSC", "CZSC-phloem-xylem", "SC-ii", "SC-oi", "SC-SUS")))

dotplot_data$time <- factor(dotplot_data$time, levels = timepoints)
dotplot_data$features.plot <- factor(dotplot_data$features.plot, levels = markers)

# ----------------------------------------------------------------------------
# 5. Generate the faceted dotplot
# ----------------------------------------------------------------------------
p <- ggplot(dotplot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
  geom_point(alpha = 0.6) +
  theme(text = element_text(size = 13),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1, size = 11),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "gray40")) +
  scale_size_area(max_size = 5) +
  labs(x = "Features", y = "Cell Clusters", size = "Percent Expressed", color = "Average Expression") +
  scale_color_gradient2(low = "blue", high = "red", mid = "white",
                        midpoint = mean(dotplot_data$avg.exp.scaled, na.rm = TRUE),
                        limit = range(dotplot_data$avg.exp.scaled, na.rm = TRUE),
                        na.value = "grey50") +
  facet_grid(celltype ~ time, scales = "free", space = "free") +
  coord_flip()

# Save the plot
ggsave("../plots/CK_markers_dotplot_top3.svg", plot = p, height = 7, width = 16)
