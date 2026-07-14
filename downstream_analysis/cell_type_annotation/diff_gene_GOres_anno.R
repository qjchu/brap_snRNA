###############################################################################
# Script: GO enrichment analysis and visualisation for cluster markers
# Steps:
#   1. Load marker genes and filter by adjusted p-value and log2FC.
#   2. Run GO enrichment (all ontologies) for each cluster.
#   3. Read combined results and create bubble plots for BP, CC, MF separately,
#      showing top 2 terms per cluster.
###############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(sva)
library(ggsci)
library(ggrepel)
library(tidyr)
library(stringr)

sessionInfo()
getwd()

if(TRUE){
  
  # Load clusterProfiler and the Brassica rapa annotation package
  library(clusterProfiler)
  library(org.Brapa.eg.db)
  
  # Read marker gene table (output from FindAllMarkers)
  marker = read.table("TT9D_res1_harmony_cluster_markers.txt", sep = "\t", header = TRUE)
  
  # Filter markers: adjusted p-value < 0.01 and absolute log2FC > 1
  marker = marker[marker$p_val_adj < 0.01,]
  marker = marker[abs(marker$avg_log2FC > 1) == 1,]
  table(marker$cluster)   # Check how many markers per cluster
  
  # ------------------- GO enrichment -------------------
  if (TRUE) {
    # Loop over clusters 1 to 22 (adjust if needed)
    for (cluster in seq(1,22)) {
      # Extract gene IDs (rownames) for the current cluster
      data1 = rownames(marker[marker$cluster == cluster,])
      
      # Run GO enrichment for all ontologies (BP, CC, MF)
      ego_ALL <- enrichGO(gene = data1, OrgDb = org.Brapa.eg.db, ont = "ALL", 
                          keyType = "GID", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      
      # If results exist, add cluster identifier and append to file
      if(!is.null(ego_ALL) && length(ego_ALL@result$Count) > 0){
        ego_ALL@result$cluster = cluster
        write.table(as.data.frame(ego_ALL@result), 
                    "TT9D_cluster_marker_GOres.txt", 
                    quote = FALSE, row.names = TRUE, col.names = TRUE, 
                    sep="\t", append = TRUE)
      }
    }
    
    # ------------------- Visualisation -------------------
    if(TRUE){
      # Read the combined GO results
      data_plot = read.table("TT9D_cluster_marker_GOres.txt", header = FALSE, sep = '\t')
      colnames(data_plot) = c('ID1', 'ONTOLOGY', 'ID2', 'Description', 
                              'GeneRatio', 'BgRatio', 'RichFactor', 
                              'FoldEnrichment', 'zScore', 'pvalue', 
                              'p.adjust', 'qvalue', 'geneID', 'Count', 'cluster')
      
      # ---- BP (Biological Process) ----
      if(TRUE){
        data_plot_BP = data_plot[data_plot$ONTOLOGY == 'BP',]
        # Take the top 2 terms per cluster (by default first 2 after arrange)
        data_plot_BP <- data_plot_BP %>% arrange(cluster) %>% group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        # Function to compute GeneRatio as numerator/denominator
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_BP$GeneRatio <- sapply(data_plot_BP$GeneRatio, divide_strings)
        data_plot_BP$Description <- str_wrap(data_plot_BP$Description, width = 50)
        
        # Bubble plot: GeneRatio vs Description, size = Count, colour = p.adjust
        ggplot(data_plot_BP, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6) +  
          scale_size_area(max_size = 10) + 
          labs(title = "TT9D", x = "GeneRatio", y = "BP", size = "Gene Count", color = "Adjust P-value") +
          scale_color_gradient(low = "blue", high = "red") +
          geom_text(aes(label = cluster), vjust = -0.5, hjust = 0.5, size = 3) +
          theme(text =  element_text(size = 15),
                panel.grid.major = element_line(colour = 'gray'), 
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = 'darkgray'), 
                panel.background = element_rect(fill = 'white'),
                plot.background = element_rect(fill="white"))
      }
      
      # ---- CC (Cellular Component) ----
      if(TRUE){
        data_plot_CC = data_plot[data_plot$ONTOLOGY == 'CC',]
        data_plot_CC <- data_plot_CC %>% arrange(cluster) %>% group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_CC$GeneRatio <- sapply(data_plot_CC$GeneRatio, divide_strings)
        data_plot_CC$Description <- str_wrap(data_plot_CC$Description, width = 50)
        
        ggplot(data_plot_CC, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6) +
          scale_size_area(max_size = 10) +
          labs(title = "TT9D", x = "GeneRatio", y = "CC", size = "Gene Count", color = "Adjust P-value") +
          scale_color_gradient(low = "blue", high = "red") +
          geom_text(aes(label = cluster), vjust = -0.5, hjust = 0.5, size = 3) +
          theme(text =  element_text(size = 15),
                panel.grid.major = element_line(colour = 'gray'), 
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = 'darkgray'), 
                panel.background = element_rect(fill = 'white'),
                plot.background = element_rect(fill="white"))
      }
      
      # ---- MF (Molecular Function) ----
      if(TRUE){
        data_plot_MF = data_plot[data_plot$ONTOLOGY == 'MF',]
        data_plot_MF <- data_plot_MF %>% arrange(cluster) %>% group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_MF$GeneRatio <- sapply(data_plot_MF$GeneRatio, divide_strings)
        data_plot_MF$Description <- str_wrap(data_plot_MF$Description, width = 50)
        
        ggplot(data_plot_MF, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6) +  
          scale_size_area(max_size = 10) + 
          labs(title = "TT9D", x = "GeneRatio", y = "MF", size = "Gene Count", color = "Adjust P-value") +
          scale_color_gradient(low = "blue", high = "red") +
          geom_text(aes(label = cluster), vjust = -0.5, hjust = 0.5, size = 3) +
          theme(text =  element_text(size = 15),
                panel.grid.major = element_line(colour = 'gray'), 
                panel.grid.minor = element_blank(), 
                axis.line = element_line(colour = 'darkgray'), 
                panel.background = element_rect(fill = 'white'),
                plot.background = element_rect(fill="white"))
      }
    }
  }
}
