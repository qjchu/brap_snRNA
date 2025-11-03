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
setwd("D:/CHU THINKBOOK/A03 白菜修改/scripts/")

# umap plot
if(TRUE){
  ## CK7D
  metadata = read.table("D:/CHU HP/A06 白菜单细胞/CK_batch3_10samples/CK7D/CK7D_metadata_res1.txt", header = TRUE, sep = "\t")
  table(metadata$harmony_clusters_res1)
  table(metadata$orig.ident)
  
  # cluster umap
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = harmony_clusters_res1)) + 
    geom_point(size = 0.9, shape = 16, alpha = 0.8) + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line =   element_blank(), 
          axis.ticks =  element_blank(), 
          axis.title =  element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'transparent'),
          plot.background=element_rect(fill = "transparent"),
          legend.position = "none") 
  ggsave('../revised_data/CK7D_res1_harmony_cluster.svg', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  
  ## TT7D
  metadata = read.table("../revised_data/TT7D/TT7D_combined_metadata_res1_5.txt", header = TRUE, sep = "\t")
  table(metadata$harmony_clusters_res1_5)
  table(metadata$orig.ident)
  
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = harmony_clusters_res1_5)) + 
    geom_point(size = 0.9, shape = 16, alpha = 0.8) + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line =   element_blank(), 
          axis.ticks =  element_blank(), 
          axis.title =  element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'transparent'),
          plot.background=element_rect(fill = "transparent"),
          legend.position = "none") 
  ggsave('../revised_data/TT7D/TT7D_combined_res1_5_harmony_cluster_nolabel.svg', plot = last_plot(), height = 8, width = 8, bg = "transparent")
} 