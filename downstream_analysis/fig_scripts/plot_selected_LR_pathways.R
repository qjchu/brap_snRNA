###############################################################################
# Purpose: Visualise selected ligand-receptor interactions across multiple
#          timepoints (CK and TT) using CellChat output.
#          For each interaction, show the source cell type, target cell type,
#          communication probability (point size), and significance (colour).
#          Summarises data from multiple CellChat edge list files.
# Input:  CellChat net_lr.txt files (containing all inferred LR interactions)
#         for each timepoint (CK3D, CK5D, ..., TT11D).
# Output: Faceted dot plot with time on x-axis, source->target on y-axis,
#         size = probability, colour = p-value.
###############################################################################

library(ggplot2)
library(ggsci)
library(dplyr)

# Set working directory (adjust to your path)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# Helper function: process a single CellChat edge file
# - Replaces numeric cluster IDs with cell type names (using recode)
# - Creates simplified source/target labels for plotting (e.g., "SC", "Endosperm")
# - Filters to keep only selected interactions (here a specific set)
# ----------------------------------------------------------------------------
process_cellchat_file <- function(file_path, time_label, source_recode, target_recode) {
  data <- read.table(file_path, header = TRUE, sep = ",")
  
  # Convert numeric source/target to cell type names
  data$source2 <- recode(data$source, !!!source_recode)
  data$target2 <- recode(data$target, !!!target_recode)
  
  # Extract main cell type (remove parenthetical details)
  data$source3 <- sub(" .*", "", data$source2)
  data$target3 <- sub(" .*", "", data$target2)
  
  # Map specific endosperm-related types to "Endosperm"
  data$source4 <- data$source3
  data$target4 <- data$target3
  data$source4 <- gsub(".*PEN.*|.*CZE.*|.*MCE.*", "Endosperm", data$source4)
  data$target4 <- gsub(".*PEN.*|.*CZE.*|.*MCE.*", "Endosperm", data$target4)
  
  # Keep only interactions of interest (selected ligand-receptor pairs)
  selected_LR <- c("ALE1->GSO2", "RALFL33->AGP16", "RALFL33->AGP20",
                   "LP1->CAM1", "CYSB->MSS1")
  plot_data <- data[data$interaction_name_2 %in% selected_LR, ]
  plot_data$time <- time_label
  return(plot_data)
}

# ----------------------------------------------------------------------------
# Define recode mappings for each timepoint (manually curated)
# ----------------------------------------------------------------------------
# CK3D (20 clusters)
ck3d_map <- setNames(c("SC-oi (c1)","SC-dividing-cell (c2)","EP-dividing-cell (c3)",
                       "SC-dividing-cell (c4)","CZSC (c5)","SC (c6)","SC-ii (c7)",
                       "SC-dividing-cell (c8)","SC (c9)","SC-endosperm (c10)",
                       "SC-ii (c11)","SC-oi (c12)","SC (c13)","SC-endosperm (c14)",
                     "SC (c15)","SC-SUS (c16)","SC (c17)","PEN/MCE (c18)",
                     "CZSC-phloem-xylem (c19)","CZE (c20)"), as.character(1:20))

# CK5D (18 clusters)
ck5d_map <- setNames(c("SC-oi (c1)","SC-endosperm (c2)","PEN (c3)",
                       "SC-dividing-cell (c4)","SC-oi (c5)","CZSC (c6)",
                       "SC (c7)","SC-ii (c8)","SC (c9)","SC-ii (c10)",
                       "EP-dividing-cell (c11)","SC-endosperm (c12)","SC (c13)",
                       "CZE (c14)","SC-SUS (c15)","MCE (c16)",
                       "CZSC-phloem-xylem (c17)","SC (c18)"), as.character(1:18))

# CK7D (25 clusters) – abbreviated for brevity; full list as in original
ck7d_map <- setNames(c("SC-oi (c1)","SC (c2)","PEN (c3)","SC-endosperm (c4)",
                       "SC-ii (c5)","SC-endosperm (c6)","SC-dividing-cell (c7)",
                       "CZSC (c8)","CZSC (c9)","SC (c10)","EP-dividing-cell (c11)",
                       "SC (c12)","SC (c13)","SC-ii (c14)","SC-endosperm (c15)",
                       "CZE (c16)","EP-dividing-cell (c17)","SC-oi (c18)","MCE (c19)",
                       "SC (c20)","SC-SUS (c21)","CZSC-phloem-xylem (c22)","SC (c23)",
                       "CZSC (c24)","EP (c25)"), as.character(1:25))

# CK9D (23 clusters)
ck9d_map <- setNames(c("SC-endosperm (c1)","CZE (c2)","SC (c3)","CZSC (c4)",
                       "SC-ii (c5)","SC-endosperm (c6)","SC (c7)","SC (c8)",
                       "PEN (c9)","EP-dividing-cell (c10)","PEN (c11)","SC-oi (c12)",
                       "SC-oi (c13)","CZSC (c14)","PEN (c15)","PEN-around-EP (c16)",
                       "SC (c17)","SC (c18)","MCE-PEN (c19)","PEN (c20)","EP (c21)",
                       "CZSC-phloem-xylem (c22)","CZE (c23)"), as.character(1:23))

# CK11D (27 clusters)
ck11d_map <- setNames(c("PEN (c1)","PEN (c2)","SC (c3)","EP (c4)","SC-ii (c5)",
                        "PEN (c6)","EP (c7)","EP (c8)","CZE (c9)","SC (c10)",
                        "SC-PEN (c11)","CZE (c12)","EP-dividing-cell (c13)",
                        "SC-oi (c14)","EP (c15)","PEN (c16)","CZSC (c17)",
                        "CZSC-phloem-xylem (c18)","CZE (c19)","CZSC (c20)",
                        "EP (c21)","PEN (c22)","CZE (c23)","EP (c24)","CZSC (c25)",
                        "SC-ii (c26)","SC (c27)"), as.character(1:27))

# TT5D (24 clusters)
tt5d_map <- setNames(c("SC (c1)","SC-endosperm (c2)","SC-oi (c3)","SC-ii (c4)",
                       "SC-dividing-cell (c5)","SC-endosperm (c6)","SC (c7)",
                       "CZSC (c8)","SC (c9)","CZSC (c10)","SC-dividing-cell (c11)",
                       "SC (c12)","SC-ii (c13)","CZSC (c14)","EP-dividing-cell (c15)",
                       "SC-SUS (c16)","CZSC (c17)","Endosperm-PCD (c18)",
                       "SC-dividing-cell (c19)","SC-dividing-cell (c20)","SC (c21)",
                       "SC (c22)","CZSC-phloem-xylem (c23)","Endosperm-active (c24)"), as.character(1:24))

# TT7D (25 clusters)
tt7d_map <- setNames(c("SC-oi (c1)","SC-ii (c2)","SC (c3)","SC-endosperm (c4)",
                       "SC (c5)","SC-ii (c6)","SC (c7)","SC (c8)","SC-dividing-cell (c9)",
                       "SC-oi (c10)","Endosperm-PCD (c11)","SC (c12)","SC (c13)",
                       "SC-ii (c14)","SC (c15)","CZSC (c16)","CZSC (c17)","SC (c18)",
                       "CZSC (c19)","SC-ii (c20)","SC-oi (c21)","EP-dividing-cell (c22)",
                       "CZSC-phloem-xylem (c23)","SC-SUS (c24)","Endosperm-active (c25)"), as.character(1:25))

# TT8D (24 clusters)
tt8d_map <- setNames(c("SC (c1)","CZSC (c2)","SC (c3)","SC-ii (c4)","SC-oi (c5)",
                       "SC-endosperm (c6)","SC-dividing-cell (c7)","SC (c8)",
                       "CZSC (c9)","SC-oi (c10)","SC (c11)","CZSC (c12)","SC-ii (c13)",
                       "SC (c14)","Endosperm-PCD (c15)","SC (c16)","EP-dividing-cell (c17)",
                       "SC (c18)","Endosperm-PCD (c19)","SC (c20)","CZSC-phloem-xylem (c21)",
                       "Endosperm-active (c22)","SC-SUS (c23)","Endosperm-active (c24)"), as.character(1:24))

# TT9D (28 clusters)
tt9d_map <- setNames(c("SC (c1)","SC (c2)","SC (c3)","SC-endosperm (c4)","CZSC (c5)",
                       "SC (c6)","SC-ii (c7)","SC-oi (c8)","SC-oi (c9)","SC-ii (c10)",
                       "SC (c11)","CZSC (c12)","CZSC (c13)","SC (c14)","CZSC (c15)",
                       "SC (c16)","SC (c17)","SC-ii (c18)","SC-dividing-cell (c19)",
                       "Endosperm-PCD (c20)","CZSC (c21)","CZSC-phloem-xylem (c22)",
                       "SC (c23)","EP-dividing-cell (c24)","SC (c25)","SC-SUS (c26)",
                       "SC (c27)","Endosperm-active (c28)"), as.character(1:28))

# TT11D (17 clusters)
tt11d_map <- setNames(c("CZSC (c1)","SC (c2)","SC (c3)","SC (c4)","SC-oi (c5)",
                        "SC (c6)","SC (c7)","SC (c8)","CZSC (c9)","SC (c10)",
                        "CZSC (c11)","SC (c12)","CZSC-phloem-xylem (c13)","SC (c14)",
                        "SC (c15)","SC (c16)","Endosperm (c17)"), as.character(1:17))

# ----------------------------------------------------------------------------
# Process all files (adjust paths to your actual data location)
# ----------------------------------------------------------------------------
base_dir <- "../path/to/cellchat_files/"

plot_list <- list()

# CK
plot_list[[1]] <- process_cellchat_file(paste0(base_dir, "CK3D_net_lr.txt"), "CK3D", ck3d_map, ck3d_map)
plot_list[[2]] <- process_cellchat_file(paste0(base_dir, "CK5D_net_lr.txt"), "CK5D", ck5d_map, ck5d_map)
plot_list[[3]] <- process_cellchat_file(paste0(base_dir, "CK7D_net_lr.txt"), "CK7D", ck7d_map, ck7d_map)
plot_list[[4]] <- process_cellchat_file(paste0(base_dir, "CK9D_net_lr.txt"), "CK9D", ck9d_map, ck9d_map)
plot_list[[5]] <- process_cellchat_file(paste0(base_dir, "CK11D_net_lr.txt"), "CK11D", ck11d_map, ck11d_map)

# TT
plot_list[[6]] <- process_cellchat_file(paste0(base_dir, "TT5D_net_lr.txt"), "TT5D", tt5d_map, tt5d_map)
plot_list[[7]] <- process_cellchat_file(paste0(base_dir, "TT7D_net_lr.txt"), "TT7D", tt7d_map, tt7d_map)
plot_list[[8]] <- process_cellchat_file(paste0(base_dir, "TT8D_net_lr.txt"), "TT8D", tt8d_map, tt8d_map)
plot_list[[9]] <- process_cellchat_file(paste0(base_dir, "TT9D_net_lr.txt"), "TT9D", tt9d_map, tt9d_map)
plot_list[[10]] <- process_cellchat_file(paste0(base_dir, "TT11D_net_lr.txt"), "TT11D", tt11d_map, tt11d_map)

plot_data <- bind_rows(plot_list)

# Add treatment (CK/TT) and source->target label
plot_data$treat <- substr(plot_data$time, 1, 2)
plot_data$sourceTarget <- paste(plot_data$source4, "->", plot_data$target4, sep = "")

# Reorder time factor
plot_data$time <- factor(plot_data$time,
                         levels = c("CK3D","CK5D","CK7D","CK9D","CK11D",
                                    "TT5D","TT7D","TT8D","TT9D","TT11D"))

# Optionally rename interaction labels (replace gene IDs with symbols)
plot_data$interaction_name_3 <- plot_data$interaction_name_2
plot_data$interaction_name_3 <- gsub("AT1G47710","Bch10G007420", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AT3G16920","Bch08G039110", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AT3G52370","Bch04G008210", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AtCDC48B","CDC48B", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AtCDC48C","CDC48C", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AT1G66250","Bch02G019130", plot_data$interaction_name_3)
plot_data$interaction_name_3 <- gsub("AT3G13560","Bch05G039790", plot_data$interaction_name_3)

# ----------------------------------------------------------------------------
# Generate plot: points show probability (size) and p-value (colour)
# Facet by interaction (ligand-receptor pair)
# ----------------------------------------------------------------------------
p <- ggplot(plot_data, aes(x = time, y = sourceTarget, size = prob, color = pval)) +
  geom_point(alpha = 0.6) +
  scale_size_continuous(range = c(2, 10)) +
  scale_color_gradient(low = "#228B22", high = "#DAA520") +
  facet_grid(~ interaction_name_3, scales = "free", space = "free") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        text = element_text(size = 15),
        axis.line = element_line(color = "darkgray"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"))

# Save
ggsave("../plots/cellchat_pathway.svg", p, height = 16, width = 12)
ggsave("../plots/cellchat_pathway.png", p, height = 10, width = 9, bg = "transparent")
