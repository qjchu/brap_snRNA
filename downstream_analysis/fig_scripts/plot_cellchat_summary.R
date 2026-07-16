###############################################################################
# Purpose: Summarise CellChat output for multiple timepoints (CK and TT)
#          and generate scatter plots showing source/target/self interaction
#          counts per cluster. Each point represents a cluster, with:
#            - X-axis: number of outgoing signals (source)
#            - Y-axis: number of incoming signals (target)
#            - Point size: number of self-signals (source=target)
#            - Colour: major cell type group (EP, SC, Endosperm)
#          Facets by timepoint are used to compare conditions.
# Input:  CellChat edge list files (net_lr.txt) for each sample.
# Output: SVG/PNG plots for CK and TT separately.
###############################################################################

library(ggplot2)
library(ggsci)

# Set working directory (adjust to your path)
setwd("path/to/your/analysis_directory")   # e.g., "scripts/"

# ----------------------------------------------------------------------------
# Helper function to process one CellChat edge list file
# Input:  file path, time label, cluster labels and types vectors
# Output: data frame with columns: aa_self, aa_source, aa_target, time, cluster, celltype, type
# ----------------------------------------------------------------------------
process_cellchat_file <- function(file_path, time_label, cluster_labels, cluster_types) {
  data <- read.table(file_path, header = TRUE, sep = ",")
  n_clusters <- length(cluster_labels)
  
  aa_self <- numeric(n_clusters)
  aa_source <- numeric(n_clusters)
  aa_target <- numeric(n_clusters)
  
  for (i in 1:n_clusters) {
    aa_self[i] <- length(data[data$target == i & data$source == i, 1])
    aa_source[i] <- length(data[data$source == i, 1]) - aa_self[i]
    aa_target[i] <- length(data[data$target == i, 1]) - aa_self[i]
  }
  
  df <- data.frame(aa_self, aa_source, aa_target)
  df$time <- time_label
  df$cluster <- paste0(time_label, "_", 1:n_clusters)
  df$celltype <- cluster_labels
  df$type <- cluster_types
  return(df)
}

# ----------------------------------------------------------------------------
# CK timepoints
# ----------------------------------------------------------------------------
# Define cluster labels and types for each timepoint (manually curated)
# The labels and types must match the order of clusters in the CellChat output.
# Paths to the CellChat output files (placeholders – adjust accordingly)

# CK3D
ck3d_labels <- c("SC-oi (c1)","SC-dividing-cell (c2)","EP-dividing-cell (c3)",
                 "SC-dividing-cell (c4)","CZSC (c5)","SC (c6)","SC-ii (c7)",
                 "SC-dividing-cell (c8)","SC (c9)","SC-endosperm (c10)",
                 "SC-ii (c11)","SC-oi (c12)","SC (c13)","SC-endosperm (c14)",
                 "SC (c15)","SC-SUS (c16)","SC (c17)","PEN/MCE (c18)",
                 "CZSC-phloem-xylem (c19)","CZE (c20)")
ck3d_types <- c("SC","SC","EP","SC","SC","SC","SC","SC","SC","SC",
                "SC","SC","SC","SC","SC","SC","SC","Endosperm","SC","Endosperm")

plot_data_CK <- process_cellchat_file("../path/to/cellchat/CK3D_net_lr.txt",
                                      "CK3D", ck3d_labels, ck3d_types)

# CK5D
ck5d_labels <- c("SC-oi (c1)", "SC-endosperm (c2)","PEN (c3)","SC-dividing-cell (c4)",
                 "SC-oi (c5)", "CZSC (c6)","SC (c7)", "SC-ii (c8)", "SC (c9)",
                 "SC-ii (c10)", "EP-dividing-cell (c11)", "SC-endosperm (c12)",
                 "SC (c13)", "CZE (c14)", "SC-SUS (c15)", "MCE (c16)",
                 "CZSC-phloem-xylem (c17)", "SC (c18)")
ck5d_types <- c("SC","SC","Endosperm","SC","SC","SC","SC","SC","SC",
                "SC","EP","SC","SC","Endosperm","SC","Endosperm","SC","SC")
plot_data_CK <- rbind(plot_data_CK,
                      process_cellchat_file("../path/to/cellchat/CK5D_net_lr.txt",
                                            "CK5D", ck5d_labels, ck5d_types))

# CK7D
ck7d_labels <- c("SC-oi (c1)", "SC (c2)", "PEN (c3)", "SC-endosperm (c4)",
                 "SC-ii (c5)","SC-endosperm (c6)", "SC-dividing-cell (c7)",
                 "CZSC (c8)", "CZSC (c9)", "SC (c10)", "EP-dividing-cell (c11)",
                 "SC (c12)", "SC (c13)", "SC-ii (c14)", "SC-endosperm (c15)",
                 "CZE (c16)", "EP-dividing-cell (c17)", "SC-oi (c18)", "MCE (c19)",
                 "SC (c20)", "SC-SUS (c21)","CZSC-phloem-xylem (c22)", "SC (c23)",
                 "CZSC (c24)", "EP (c25)")
ck7d_types <- c("SC","SC","Endosperm","SC","SC","SC","SC","SC","SC","SC",
                "EP","SC","SC","SC","SC","Endosperm","EP","SC","Endosperm",
                "SC","SC","SC","SC","SC","EP")
plot_data_CK <- rbind(plot_data_CK,
                      process_cellchat_file("../path/to/cellchat/CK7D_net_lr.txt",
                                            "CK7D", ck7d_labels, ck7d_types))

# CK9D
ck9d_labels <- c("SC-endosperm (c1)", "CZE (c2)", "SC (c3)", "CZSC (c4)",
                 "SC-ii (c5)","SC-endosperm (c6)", "SC (c7)", "SC (c8)",
                 "PEN (c9)", "EP-dividing-cell (c10)", "PEN (c11)", "SC-oi (c12)",
                 "SC-oi (c13)", "CZSC (c14)", "PEN (c15)", "PEN-around-EP (c16)",
                 "SC (c17)", "SC (c18)", "MCE-PEN (c19)", "PEN (c20)", "EP (c21)",
                 "CZSC-phloem-xylem (c22)", "CZE (c23)")
ck9d_types <- c("SC","Endosperm","SC","SC","SC","SC","SC","SC","Endosperm",
                "EP","Endosperm","SC","SC","SC","Endosperm","Endosperm","SC","SC",
                "Endosperm","Endosperm","EP","SC","Endosperm")
plot_data_CK <- rbind(plot_data_CK,
                      process_cellchat_file("../path/to/cellchat/CK9D_net_lr.txt",
                                            "CK9D", ck9d_labels, ck9d_types))

# CK11D
ck11d_labels <- c("PEN (c1)", "PEN (c2)", "SC (c3)", "EP (c4)","SC-ii (c5)",
                  "PEN (c6)", "EP (c7)","EP (c8)", "CZE (c9)", "SC (c10)",
                  "SC-PEN (c11)", "CZE (c12)", "EP-dividing-cell (c13)",
                  "SC-oi (c14)", "EP (c15)","PEN (c16)", "CZSC (c17)",
                  "CZSC-phloem-xylem (c18)", "CZE (c19)", "CZSC (c20)",
                  "EP (c21)", "PEN (c22)", "CZE (c23)","EP (c24)", "CZSC (c25)",
                  "SC-ii (c26)", "SC (c27)")
ck11d_types <- c("Endosperm","Endosperm","SC","EP","SC","Endosperm","EP","EP",
                 "Endosperm","SC","SC","Endosperm","EP","SC","EP","Endosperm",
                 "SC","SC","Endosperm","SC","EP","Endosperm","Endosperm","EP",
                 "SC","SC","SC")
plot_data_CK <- rbind(plot_data_CK,
                      process_cellchat_file("../path/to/cellchat/CK11D_net_lr.txt",
                                            "CK11D", ck11d_labels, ck11d_types))

# Add treatment label
plot_data_CK$treat <- "CK"

# ----------------------------------------------------------------------------
# TT timepoints (mutant) – similar structure
# ----------------------------------------------------------------------------
# TT5D
tt5d_labels <- c("SC (c1)", "SC-endosperm (c2)", "SC-oi (c3)", "SC-ii (c4)",
                 "SC-dividing-cell (c5)", "SC-endosperm (c6)", "SC (c7)",
                 "CZSC (c8)", "SC (c9)", "CZSC (c10)", "SC-dividing-cell (c11)",
                 "SC (c12)", "SC-ii (c13)", "CZSC (c14)", "EP-dividing-cell (c15)",
                 "SC-SUS (c16)", "CZSC (c17)", "Endosperm-PCD (c18)",
                 "SC-dividing-cell (c19)", "SC-dividing-cell (c20)", "SC (c21)",
                 "SC (c22)", "CZSC-phloem-xylem (c23)", "Endosperm-active (c24)")
tt5d_types <- c(rep("SC",14),"EP","SC","SC","Endosperm",rep("SC",5),"Endosperm")
plot_data_TT <- process_cellchat_file("../path/to/cellchat/TT5D_net_lr.txt",
                                      "TT5D", tt5d_labels, tt5d_types)

# TT7D
tt7d_labels <- c("SC-oi (c1)", "SC-ii (c2)", "SC (c3)", "SC-endosperm (c4)",
                 "SC (c5)", "SC-ii (c6)", "SC (c7)", "SC (c8)",
                 "SC-dividing-cell (c9)", "SC-oi (c10)", "Endosperm-PCD (c11)",
                 "SC (c12)", "SC (c13)", "SC-ii (c14)", "SC (c15)", "CZSC (c16)",
                 "CZSC (c17)", "SC (c18)", "CZSC (c19)", "SC-ii (c20)",
                 "SC-oi (c21)", "EP-dividing-cell (c22)", "CZSC-phloem-xylem (c23)",
                 "SC-SUS (c24)", "Endosperm-active (c25)")
tt7d_types <- c(rep("SC",10),"Endosperm",rep("SC",10),"EP","SC","SC","Endosperm")
plot_data_TT <- rbind(plot_data_TT,
                      process_cellchat_file("../path/to/cellchat/TT7D_net_lr.txt",
                                            "TT7D", tt7d_labels, tt7d_types))

# TT8D
tt8d_labels <- c("SC (c1)", "CZSC (c2)", "SC (c3)", "SC-ii (c4)", "SC-oi (c5)",
                 "SC-endosperm (c6)", "SC-dividing-cell (c7)", "SC (c8)",
                 "CZSC (c9)", "SC-oi (c10)", "SC (c11)", "CZSC (c12)",
                 "SC-ii (c13)", "SC (c14)", "Endosperm-PCD (c15)", "SC (c16)",
                 "EP-dividing-cell (c17)", "SC (c18)", "Endosperm-PCD (c19)",
                 "SC (c20)", "CZSC-phloem-xylem (c21)", "Endosperm-active (c22)",
                 "SC-SUS (c23)", "Endosperm-active (c24)")
tt8d_types <- c(rep("SC",14),"Endosperm","SC","EP","SC","Endosperm","SC","SC",
                "Endosperm","SC","Endosperm")
plot_data_TT <- rbind(plot_data_TT,
                      process_cellchat_file("../path/to/cellchat/TT8D_net_lr.txt",
                                            "TT8D", tt8d_labels, tt8d_types))

# TT9D
tt9d_labels <- c("SC (c1)", "SC (c2)", "SC (c3)", "SC-endosperm (c4)", "CZSC (c5)",
                 "SC (c6)", "SC-ii (c7)", "SC-oi (c8)", "SC-oi (c9)", "SC-ii (c10)",
                 "SC (c11)", "CZSC (c12)", "CZSC (c13)", "SC (c14)", "CZSC (c15)",
                 "SC (c16)", "SC (c17)", "SC-ii (c18)", "SC-dividing-cell (c19)",
                 "Endosperm-PCD (c20)", "CZSC (c21)", "CZSC-phloem-xylem (c22)",
                 "SC (c23)", "EP-dividing-cell (c24)", "SC (c25)", "SC-SUS (c26)",
                 "SC (c27)", "Endosperm-active (c28)")
tt9d_types <- c(rep("SC",19),"Endosperm","SC","SC","SC","EP","SC","SC","SC","Endosperm")
plot_data_TT <- rbind(plot_data_TT,
                      process_cellchat_file("../path/to/cellchat/TT9D_net_lr.txt",
                                            "TT9D", tt9d_labels, tt9d_types))

# TT11D
tt11d_labels <- c("CZSC (c1)", "SC (c2)", "SC (c3)", "SC (c4)", "SC-oi (c5)",
                  "SC (c6)", "SC (c7)", "SC (c8)", "CZSC (c9)", "SC (c10)",
                  "CZSC (c11)", "SC (c12)", "CZSC-phloem-xylem (c13)", "SC (c14)",
                  "SC (c15)", "SC (c16)", "Endosperm (c17)")
tt11d_types <- c(rep("SC",16), "Endosperm")
plot_data_TT <- rbind(plot_data_TT,
                      process_cellchat_file("../path/to/cellchat/TT11D_net_lr.txt",
                                            "TT11D", tt11d_labels, tt11d_types))

plot_data_TT$treat <- "TT"

# ----------------------------------------------------------------------------
# Generate plots
# ----------------------------------------------------------------------------
# CK plot
plot_data_CK$time <- factor(plot_data_CK$time, levels = c("CK3D", "CK5D", "CK7D", "CK9D", "CK11D"))
p_ck <- ggplot(plot_data_CK, aes(x = aa_source, y = aa_target, size = aa_self, color = type)) +
  geom_point() +
  geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 1.8) +
  scale_color_npg("nrc") +
  scale_size_continuous(range = c(2, 15)) +
  xlim(0, 350) + ylim(0, 350) +
  facet_grid(~ time, scales = "free", space = "free") +
  theme(text = element_text(size = 15),
        axis.line = element_line(color = "darkgray"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "right")
ggsave("../plots/CK_cellchat.svg", p_ck, height = 4.3, width = 20)
ggsave("../plots/CK_cellchat.png", p_ck, height = 4.3, width = 20, bg = "transparent")

# TT plot
plot_data_TT$time <- factor(plot_data_TT$time, levels = c("TT5D", "TT7D", "TT8D", "TT9D", "TT11D"))
p_tt <- ggplot(plot_data_TT, aes(x = aa_source, y = aa_target, size = aa_self, color = type)) +
  geom_point() +
  geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 1.8) +
  scale_color_npg("nrc") +
  scale_size_continuous(range = c(2, 12)) +
  xlim(0, 350) + ylim(0, 350) +
  facet_grid(~ time, scales = "free", space = "free") +
  theme(text = element_text(size = 15),
        axis.line = element_line(color = "darkgray"),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        legend.position = "right")
ggsave("../plots/TT_cellchat.svg", p_tt, height = 4.3, width = 20)
ggsave("../plots/TT_cellchat.png", p_tt, height = 4.3, width = 20, bg = "transparent")
