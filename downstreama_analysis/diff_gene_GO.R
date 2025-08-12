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
  
  library(clusterProfiler)
  library(org.Brapa.eg.db)
  
  marker = read.table("TT9D_res1_harmony_cluster_markers.txt", sep = "\t", header = TRUE)
  marker = marker[marker$p_val_adj < 0.01,]
  marker = marker[abs(marker$avg_log2FC > 1) == 1,]
  table(marker$cluster)
  
  # GO
  if (TRUE) {
    for (cluster in seq(1,22)) {
      data1 = rownames(marker[marker$cluster == cluster,])
      # bitr(data1, 'GID', c("EVIDENCE", "EVIDENCEALL", "GENENAME", "GO", "GOALL", "Ko", "ONTOLOGY", "ONTOLOGYALL"), org.Brapa.eg.db)
      
      ego_ALL <- enrichGO(gene = data1, OrgDb = org.Brapa.eg.db, ont = "ALL", 
                          keyType = "GID", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
      
      if(!is.null(ego_ALL) && length(ego_ALL@result$Count) > 0){
        ego_ALL@result$cluster = cluster
        write.table(as.data.frame(ego_ALL@result), "TT9D_cluster_marker_GOres.txt", quote = FALSE, row.names = TRUE, col.names = TRUE, sep="\t", append = TRUE)
      }
    }
    
    # plot
    if(TRUE){
      data_plot = read.table("TT9D_cluster_marker_GOres.txt", header = FALSE, sep = '\t')
      colnames(data_plot) = c('ID1', 'ONTOLOGY',	'ID2',	'Description', 'GeneRatio',	'BgRatio',	'RichFactor',	'FoldEnrichment',	'zScore',	'pvalue',	'p.adjust',	'qvalue',	'geneID',	'Count', 'cluster')
      
      # BP
      if(TRUE){
        data_plot_BP = data_plot[data_plot$ONTOLOGY == 'BP',]
        # data_plot_BP <- data_plot_BP[!duplicated(data_plot_BP$cluster), ]
        data_plot_BP <- data_plot_BP %>% arrange(cluster) %>%  group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_BP$GeneRatio <- sapply(data_plot_BP$GeneRatio, divide_strings)
        data_plot_BP$Description <- str_wrap(data_plot_BP$Description, width = 50)
        
        ggplot(data_plot_BP, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6,) +  
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
      
      # CC
      if(TRUE){
        data_plot_CC = data_plot[data_plot$ONTOLOGY == 'CC',]
        data_plot_CC <- data_plot_CC %>% arrange(cluster) %>%  group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_CC$GeneRatio <- sapply(data_plot_CC$GeneRatio, divide_strings)
        data_plot_CC$Description <- str_wrap(data_plot_CC$Description, width = 50)
        
        ggplot(data_plot_CC, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6,) +
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
      
      # MF
      if(TRUE){
        data_plot_MF = data_plot[data_plot$ONTOLOGY == 'MF',]
        data_plot_MF <- data_plot_MF %>% arrange(cluster) %>%  group_by(cluster) %>% slice_head(n = 2) %>% ungroup()
        
        divide_strings <- function(str) {
          parts <- strsplit(str, "/")[[1]]
          numerator <- as.numeric(parts[1])
          denominator <- as.numeric(parts[2])
          if (denominator == 0) { return(NA) } else { return(numerator / denominator) }
        }
        data_plot_MF$GeneRatio <- sapply(data_plot_MF$GeneRatio, divide_strings)
        
        data_plot_MF$Description <- str_wrap(data_plot_MF$Description, width = 50)
        
        ggplot(data_plot_MF, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
          geom_point(alpha = 0.6,) +  
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
