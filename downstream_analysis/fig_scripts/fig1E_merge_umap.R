library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(ggsci)
library(ggrepel)
library(tidyr)

getwd()
setwd("scripts/")

# time-point umap plot CK
if(TRUE){
  metadata = read.table("../CK_batch3_10samples/10samples/CK_10samples_batch3_metadata_res1.txt", header = TRUE, sep = "\t")
  table(metadata$harmony_clusters_res1)
  table(metadata$orig.ident)
  
  metadata$orig.ident2 = gsub('_.*', '', metadata$orig.ident)
  celltype_med <- metadata %>% group_by(orig.ident2) %>% summarise( umapx.harmony = mean(umapx.harmony), umapy.harmony = median(umapy.harmony))
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.45, shape = 16, alpha = 0.8) + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line =   element_blank(), 
          axis.ticks =  element_blank(), 
          axis.title =  element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background=element_rect(fill="white"),
          legend.position = "none") + 
    facet_wrap(.~orig.ident2, ncol = 5)
  # ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_origident2.svg', plot = last_plot(), height = 7, width = 32)
  # ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_origident2.png', plot = last_plot(), height = 7, width = 32, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'CK3D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#9ACD32') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  # ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK3D.svg', plot = last_plot(), height = 8, width = 8)
  ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK3D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'CK5D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#6B8E23') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  # ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK5D.svg', plot = last_plot(), height = 8, width = 8)
  ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK5D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'CK7D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#556B2F') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK7D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'CK9D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#32CD32') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK9D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'CK11D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#228B22') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../CK_batch3_10samples/10samples/10samples_batch3_res1_umap_CK11D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  if(TRUE){
    table(metadata[metadata$harmony_clusters_res1 == 1,]$orig.ident)
    
    cluster_sample_counts <- metadata %>%
      dplyr::group_by(harmony_clusters_res1, orig.ident2) %>%  
      dplyr::summarise(count = n(), .groups = 'drop') %>%  
      dplyr::ungroup()
    
    sum(cluster_sample_counts[cluster_sample_counts$orig.ident2 == 'CK11D','count'])
    cluster_sample_counts$perc = cluster_sample_counts$count/18316
    a = cluster_sample_counts[cluster_sample_counts$orig.ident2 == 'CK11D',]
    order(-a$perc)
    as.vector(a[order(-a$perc), 1])
    
    cluster_sample_long <- cluster_sample_counts %>%
      pivot_longer(cols = c(count), names_to = "metric", values_to = "value") %>%
      mutate(metric = factor(metric, levels = c("count")))
    
    cluster_sample_long$harmony_clusters_res1 = factor(cluster_sample_long$harmony_clusters_res1, levels = c(seq(1,31)))
    ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, fill = orig.ident2)) +
      geom_bar(stat = "identity", position = "fill", aes(group = orig.ident2), width = 0.6, alpha = 0.6) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(text = element_text(size = 15),
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            axis.line = element_line(colour = "black",lineend = 'round'),
            legend.position = "right") +
      scale_fill_igv("default") + 
      # scale_fill_manual(values = c("TT9D_1" = "#cff09e", "TT9D_2" = "#a8dba8", 
      #                              "TT9D_1" = "#00dffc", "TT9D_2" = "#008c9e",
      #                              "TT9D_1" = "#D499B9", "TT9D_2" = "#9055A2",
      #                              "CK9D_1" = "#fec9c9", "CK9D_2" = "#ee6e9f",
      #                              "TT9D_1" = "#ede574", "TT9D_2" = "#f9d423",
      #                              "TT9D_1" = "#F6B352", "TT9D_2" = "#e3632d", "TT9D_2" = "#44633F")) +
      labs(y = "Proportion", fill = "Samples")
    
    table(cluster_sample_long$harmony_clusters_res1)
    cluster_sample_long[cluster_sample_long$harmony_clusters_res1 == 18,]
    cluster_sample_long[cluster_sample_long$harmony_clusters_res1 == 22,]
    cluster_sample_long[cluster_sample_long$harmony_clusters_res1 == 24,]
    cluster_sample_long[cluster_sample_long$harmony_clusters_res1 == 30,]
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 18, orig.ident2 = 'CK11D', perc = 0, metric = 'count', value = 0))
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 22, orig.ident2 = 'CK3D', perc = 0, metric = 'count', value = 0))
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 24, orig.ident2 = 'CK3D', perc = 0, metric = 'count', value = 0))
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 30, orig.ident2 = 'CK11D', perc = 0, metric = 'count', value = 0))
    
    cluster_sample_long$orig.ident2 = factor(cluster_sample_long$orig.ident2, levels = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D'))
    cluster_sample_long$harmony_clusters_res1 = factor(cluster_sample_long$harmony_clusters_res1,
                                                       levels = rev(c(18,30,2,21,26,16,20,29,28,3,1,13,31,5,4,6,27,12,7,11,23,22,15,10,8,9,14,17,19,24,25)))
    ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, color = orig.ident2, fill = orig.ident2)) +
      geom_area(stat = "identity", position = "fill", alpha = .8, aes(group = orig.ident2), na.rm = FALSE) +
      # scale_fill_manual(values = c('white','white','white','white','white')) + # '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22' 
      # scale_color_manual(values = c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22')) +
      scale_fill_manual(values = c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22')) +
      scale_color_manual(values = c('white','white','white','white','white')) + 
      theme(text = element_text(size = 15),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            axis.line = element_line(colour = "black",lineend = 'round'),
            legend.position = "right") +
      coord_flip() 
    ggsave('../revised_data/CK_10samples_batch3_cell_perc.svg', plot = last_plot(), height = 7, width = 3, bg = "transparent")
    
  }
}

# time-point umap plot TT
if(TRUE){
  metadata = read.table("../revised_data/TT_8samples/TT_8samples_combined_batch3_metadata_res1.txt", header = TRUE, sep = "\t")
  table(metadata$harmony_clusters_res1)
  table(metadata$orig.ident)
  
  # orig.ident umap2
  metadata$orig.ident2 = gsub('_.*', '', metadata$orig.ident)
  celltype_med <- metadata %>% group_by(orig.ident2) %>% summarise( umapx.harmony = mean(umapx.harmony), umapy.harmony = median(umapy.harmony))
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.45, shape = 16, alpha = 0.8) + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line =   element_blank(), 
          axis.ticks =  element_blank(), 
          axis.title =  element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background=element_rect(fill="white"),
          legend.position = "none") + 
    facet_wrap(.~orig.ident2, ncol = 5)
  # ggsave('../revised_data/TT_8samples/TT_combined_res1_origident2.svg', plot = last_plot(), height = 7, width = 32)
  # ggsave('../revised_data/TT_8samples/TT_combined_res1_origident2.png', plot = last_plot(), height = 7, width = 32, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'TT5D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#BDB76B') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../revised_data/TT_8samples/TT_combined_res1_umap_TT5D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'TT7D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#FFD700') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../revised_data/TT_8samples/TT_combined_res1_umap_TT7D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'TT8D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#FFA500') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../revised_data/TT_8samples/TT_combined_res1_umap_TT8D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'TT9D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#DAA520') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../revised_data/TT_8samples/TT_combined_res1_umap_TT9D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  ggplot(metadata[metadata$orig.ident2 == 'TT11D',], aes(umapx.harmony, umapy.harmony, color = orig.ident2)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98, color = '#B8860B') + xlab('UMAP1') + ylab('UMAP2') + 
    scale_color_igv("default") + 
    theme(text =  element_text(size = 15), panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),  axis.line = element_blank(),  axis.ticks =  element_blank(), axis.title =  element_blank(),  axis.text = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent"), legend.position = "none")
  ggsave('../revised_data/TT_8samples/TT_combined_res1_umap_TT11D.png', plot = last_plot(), height = 8, width = 8, bg = "transparent")
  
  if(TRUE){
    table(metadata[metadata$harmony_clusters_res1 == 1,]$orig.ident)
    head(metadata)
    
    cluster_sample_counts <- metadata %>%
      dplyr::group_by(harmony_clusters_res1, orig.ident2) %>%  
      dplyr::summarise(count = n(), .groups = 'drop') %>%  
      dplyr::ungroup()
    
    sum(cluster_sample_counts[cluster_sample_counts$orig.ident2 == 'TT11D','count'])
    cluster_sample_counts$perc = cluster_sample_counts$count/5146
    a = cluster_sample_counts[cluster_sample_counts$orig.ident2 == 'TT11D',]
    order(-a$perc)
    as.vector(a[order(-a$perc), 1])
    
    cluster_sample_long <- cluster_sample_counts %>%
      pivot_longer(cols = c(count), names_to = "metric", values_to = "value") %>%
      mutate(metric = factor(metric, levels = c("count")))
    
    cluster_sample_long$harmony_clusters_res1 = factor(cluster_sample_long$harmony_clusters_res1,
                                                       levels = c(seq(1,27)))
    ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, fill = orig.ident2)) +
      geom_bar(stat = "identity", position = "fill", aes(group = orig.ident2), width = 0.6, alpha = 0.6) +
      scale_y_continuous(labels = scales::percent_format()) +
      theme(text = element_text(size = 15),
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            axis.line = element_line(colour = "black",lineend = 'round'),
            legend.position = "right") +
      scale_fill_igv("default") + 
      # scale_fill_manual(values = c("TT9D_1" = "#cff09e", "TT9D_2" = "#a8dba8", 
      #                              "TT9D_1" = "#00dffc", "TT9D_2" = "#008c9e",
      #                              "TT9D_1" = "#D499B9", "TT9D_2" = "#9055A2",
      #                              "CK9D_1" = "#fec9c9", "CK9D_2" = "#ee6e9f",
      #                              "TT9D_1" = "#ede574", "TT9D_2" = "#f9d423",
      #                              "TT9D_1" = "#F6B352", "TT9D_2" = "#e3632d", "TT9D_2" = "#44633F")) +
      labs(y = "Proportion", fill = "Samples")
    
    table(cluster_sample_long$harmony_clusters_res1)
    cluster_sample_long[cluster_sample_long$harmony_clusters_res1 == 27,]
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 27, orig.ident2 = 'TT9D', perc = 0, metric = 'count', value = 0))
    cluster_sample_long = rbind(cluster_sample_long, data.frame(harmony_clusters_res1 = 27, orig.ident2 = 'TT11D', perc = 0, metric = 'count', value = 0))
    
    cluster_sample_long$orig.ident2 = factor(cluster_sample_long$orig.ident2, levels = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D'))
    cluster_sample_long$harmony_clusters_res1 = factor(cluster_sample_long$harmony_clusters_res1,
                                                       levels = c(22,20,2,8,3,16,24,14,4,23,26,5,17,19,7,25,15,18,6,1,13,9,10,12,11,21,27))
    # ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, color = orig.ident2, fill = orig.ident2)) +
    #   geom_area(stat = "identity", position = "fill", alpha = .8, aes(group = orig.ident2), na.rm = FALSE) +
    #   # scale_fill_manual(values = c('white','white','white','white','white')) + # '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22' 
    #   # scale_color_manual(values = c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22')) +
    #   scale_fill_manual(values = c( '#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
    #   scale_color_manual(values = c('white','white','white','white','white')) + 
    #   theme(text = element_text(size = 15),
    #         axis.title.y = element_blank(),
    #         panel.background = element_rect(fill = 'white'),
    #         plot.background=element_rect(fill="white"),
    #         axis.line = element_line(colour = "black",lineend = 'round'),
    #         legend.position = "right") 
    # ggsave('../revised_data/TT_8samples/TT_combined_cell_perc.svg', plot = last_plot(), height = 3, width = 8.5, bg = "transparent")
    
    ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, color = orig.ident2, fill = orig.ident2)) +
      geom_area(stat = "identity", position = "fill", alpha = .8, aes(group = orig.ident2), na.rm = FALSE) +
      # scale_fill_manual(values = c('white','white','white','white','white')) + # '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22' 
      # scale_color_manual(values = c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22')) +
      scale_fill_manual(values = c('#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
      scale_color_manual(values = c('white','white','white','white','white')) + 
      theme(text = element_text(size = 15),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.background = element_rect(fill = 'white'),
            plot.background=element_rect(fill="white"),
            axis.line = element_line(colour = "black",lineend = 'round'),
            legend.position = "right") +
      coord_flip() 
    ggsave('../revised_data/TT_8samples/TT_combined_cell_perc.svg', plot = last_plot(), height = 7, width = 3, bg = "transparent")
    
  }
}

