###############################################################################
# Script: UMAP visualisation and marker gene dot plots for snRNA-seq data
#         (TT9D combined dataset)
# Steps:
#   1. Load metadata and generate UMAP plots (cluster labels & sample origin)
#   2. Calculate and plot sample proportion per cluster
#   3. Load Seurat object and generate dot plots for various gene sets:
#       - PlantscRNAdb markers (embryo, seed coat, endosperm)
#       - Region-specific markers (pg, g, h stages)
#       - Subregion-specific markers
#       - qPCR-validated markers
#       - Cell cycle markers
#       - FISH-validated markers
###############################################################################

# Load required libraries
library(dplyr)          # Data manipulation
library(Seurat)         # Single-cell analysis
library(patchwork)      # (not actively used)
library(ggplot2)        # Plotting
library(svglite)        # SVG export
library(sva)            # (not actively used)
library(ggsci)          # Scientific colour palettes (e.g., scale_color_igv)
library(ggrepel)        # Label repulsion for geom_text_repel
library(tidyr)          # Data reshaping (pivot_longer)
library(stringr)        # String manipulation

# Print session info and working directory for reproducibility
sessionInfo()
getwd()

# ------------------- Part 1: UMAP visualisation -------------------
if(TRUE){
  # Load metadata (contains cluster assignments and embedding coordinates)
  metadata = read.table("TT9D_combined_metadata_res1_Louvain.txt", header = TRUE, sep = "\t")
  table(metadata$harmony_clusters_res1)   # Check cluster distribution
  table(metadata$orig.ident)              # Check sample distribution
  
  # 1a. UMAP coloured by Harmony cluster
  metadata$harmony_clusters_res1 = factor(metadata$harmony_clusters_res1, levels = c(seq(0,27)))
  # Calculate median/mean UMAP coordinates per cluster for label placement
  celltype_med <- metadata %>% 
    group_by(harmony_clusters_res1) %>% 
    summarise(umapx.harmony = mean(umapx.harmony), umapy.harmony = median(umapy.harmony))
  
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = harmony_clusters_res1)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.98) +
    scale_color_igv("default") + 
    theme(text = element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill = "white"),
          legend.position = "none") +
    geom_label_repel(data = celltype_med, aes(label = harmony_clusters_res1), 
                     fontface = "bold", point.padding = unit(0.5, "lines"))
  ggsave('TT9D_combined_res1_Louvain_harmony_cluster.svg', plot = last_plot(), height = 8, width = 8)
  
  # 1b. UMAP coloured by original sample (orig.ident)
  celltype_med <- metadata %>% 
    group_by(orig.ident) %>% 
    summarise(umapx.harmony = mean(umapx.harmony), umapy.harmony = median(umapy.harmony))
  
  ggplot(metadata, aes(umapx.harmony, umapy.harmony, color = orig.ident)) + 
    geom_point(size = 0.5, shape = 16, alpha = 0.6) +
    scale_color_igv("default") + 
    theme(text = element_text(size = 15),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill = "white"),
          legend.position = "none") +
    geom_label_repel(data = celltype_med, aes(label = orig.ident), 
                     fontface = "bold", point.padding = unit(0.5, "lines"))
  ggsave('TT9D_combined_res1_Louvain_origident.svg', plot = last_plot(), height = 8, width = 8)
}

# ------------------- Part 2: Sample proportion per cluster -------------------
if(TRUE){
  # Count cells per cluster per sample
  cluster_sample_counts <- metadata %>%
    group_by(harmony_clusters_res1, orig.ident) %>%  
    summarise(count = n(), .groups = 'drop') %>%  
    ungroup()
  
  # Reshape for stacked bar plot
  cluster_sample_long <- cluster_sample_counts %>%
    pivot_longer(cols = c(count), names_to = "metric", values_to = "value") %>%
    mutate(metric = factor(metric, levels = c("count")))
  cluster_sample_long$harmony_clusters_res1 = factor(cluster_sample_long$harmony_clusters_res1, levels = c(seq(1,27)))
  
  # Stacked bar chart (fill = proportion)
  ggplot(cluster_sample_long, aes(x = harmony_clusters_res1, y = value, fill = orig.ident)) +
    geom_bar(stat = "identity", position = "fill", aes(group = orig.ident), width = 0.6, alpha = 0.6) +
    scale_y_continuous(labels = scales::percent_format()) +
    theme(text = element_text(size = 15),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black", lineend = 'round'),
          legend.position = "right") +
    scale_fill_igv("default") + 
    labs(y = "Proportion", fill = "Samples")
}

# ------------------- Part 3: Marker gene dot plots -------------------
if(TRUE){
  # Load the Seurat object
  sce = readRDS("TT9D_merge_res1.rds")
  group = 'harmony_clusters_res1'
  
  ## 3a. PlantscRNAdb markers (embryo, seed coat, endosperm)
  if(TRUE){
    # Embryo markers
    DotPlot(sce, group.by = group, features = c('Bch06G016650','Bch07G018330', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Seed coat markers
    DotPlot(sce, group.by = group, features = c('Bch03G020860','Bch09G045940', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # Endosperm markers
    DotPlot(sce, group.by = group, features = c('Bch03G001460','Bch02G001240', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }
  
  ## 3b. Region-specific markers (pg, g, h stages) – subdivided by tissue
  if(TRUE){
    # pgEndosperm
    DotPlot(sce, group.by = group, features = c('Bch09G015880','Bch02G029920', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # pgSeedCoat (long list, with smaller text)
    DotPlot(sce, group.by = group, features = c('Bch02G044780','Bch09G003830', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gEmbryo
    DotPlot(sce, group.by = group, features = c('Bch03G001790','Bch05G011170', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gEndosperm
    DotPlot(sce, group.by = group, features = c('Bch01G046510','Bch05G042410')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # gSeedCoat (long)
    DotPlot(sce, group.by = group, features = c('Bch01G051440','Bch05G048410', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hEndosperm
    DotPlot(sce, group.by = group, features = c('Bch06G018610','Bch07G040550', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    # hSeedCoat (long)
    DotPlot(sce, group.by = group, features = c('Bch02G004890','Bch05G048750', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional hSeedCoat (long)
    DotPlot(sce, group.by = group, features = c('Bch01G051440','Bch05G048410', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
  }
  
  ## 3c. Subregion-specific markers (pg, g, h stages with finer subdivisions)
  if(TRUE){
    # pgEP
    DotPlot(sce, group.by = group, features = c('Bch04G036100','Bch05G005530', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # pgMCE
    DotPlot(sce, group.by = group, features = c('Bch02G027190','Bch07G030070', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # pgPEN
    DotPlot(sce, group.by = group, features = c('Bch01G023990','Bch04G032380', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # pgCZE
    DotPlot(sce, group.by = group, features = c('Bch05G002030','Bch04G033420', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional pgCZE
    DotPlot(sce, group.by = group, features = c('Bch06G013850','Bch07G010460', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # pgCZSC
    DotPlot(sce, group.by = group, features = c('Bch03G027350','Bch05G006740', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional pgCZSC
    DotPlot(sce, group.by = group, features = c('Bch01G014720','Bch01G012920', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # pgSC
    DotPlot(sce, group.by = group, features = c('Bch02G004930','Bch02G004910', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional pgSC
    DotPlot(sce, group.by = group, features = c('Bch08G033990','Bch07G016810', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    
    # gEP
    DotPlot(sce, group.by = group, features = c('Bch09G042150','Bch07G013960', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gSUS
    DotPlot(sce, group.by = group, features = c('Bch02G005010','Bch05G026550', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional gSUS
    DotPlot(sce, group.by = group, features = c('Bch09G044030','Bch06G010650', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gMCE
    DotPlot(sce, group.by = group, features = c('Bch05G002470','Bch05G002480', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gPEN
    DotPlot(sce, group.by = group, features = c('Bch04G032380','Bch03G022550', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gCZE
    DotPlot(sce, group.by = group, features = c('Bch07G036750','Bch05G017230', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional gCZE
    DotPlot(sce, group.by = group, features = c('Bch03G052180','Bch08G016210', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gCZSC
    DotPlot(sce, group.by = group, features = c('Bch05G006740','Bch03G022640', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional gCZSC
    DotPlot(sce, group.by = group, features = c('Bch02G033360','Bch01G028130', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # gSC
    DotPlot(sce, group.by = group, features = c('Bch02G004930','Bch02G004910', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional gSC
    DotPlot(sce, group.by = group, features = c('Bch09G074580','Bch08G033990', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    
    # hEP
    DotPlot(sce, group.by = group, features = c('Bch03G027320','Bch08G009720', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hMCE
    DotPlot(sce, group.by = group, features = c('Bch10G028740','Bch03G006050', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hPEN
    DotPlot(sce, group.by = group, features = c('Bch07G027820','Bch04G037090', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hCZE
    DotPlot(sce, group.by = group, features = c('Bch05G002460','Bch03G023660', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional hCZE (multiple parts)
    DotPlot(sce, group.by = group, features = c('Bch04G021690','Bch02G003740', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Further hCZE
    DotPlot(sce, group.by = group, features = c('Bch08G035430','Bch06G009620', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hCZSC
    DotPlot(sce, group.by = group, features = c('Bch05G006740','Bch03G022640', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional hCZSC
    DotPlot(sce, group.by = group, features = c('Bch08G006070','Bch01G014920', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # hSC
    DotPlot(sce, group.by = group, features = c('Bch08G014230','Bch01G020570', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Additional hSC
    DotPlot(sce, group.by = group, features = c('Bch08G021800','Bch01G007660', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
    # Further hSC
    DotPlot(sce, group.by = group, features = c('Bch08G017400','Bch01G012600', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 2))
  }
  
  ## 3d. qPCR-validated markers (shorter lists)
  if(TRUE){
    # Preglobular - Peripheral Endosperm
    DotPlot(sce, group.by = group, features = c('Bch10G003150','Bch08G040130', ...)) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # pgSC
    DotPlot(sce, group.by = group, features = c('Bch01G045500','Bch05G039810')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # pgEP
    DotPlot(sce, group.by = group, features = c('Bch06G047330','Bch07G041900','Bch02G024660')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # pgCZSC
    DotPlot(sce, group.by = group, features = c('Bch01G000700','Bch06G047560','Bch09G024150')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gCZE
    DotPlot(sce, group.by = group, features = c('Bch08G017430','Bch01G027140','Bch06G020540','Bch06G028740')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gCZSC
    DotPlot(sce, group.by = group, features = c('Bch07G043360','Bch05G030960','Bch01G037420','Bch01G015740','Bch03G055910','Bch06G018290','Bch08G006070','Bch08G019980','Bch08G019990','Bch02G008290','Bch03G009060','Bch10G024700')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gMCE
    DotPlot(sce, group.by = group, features = c('Bch03G029600')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gPEN
    DotPlot(sce, group.by = group, features = c('Bch10G033410')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gSC
    DotPlot(sce, group.by = group, features = c('Bch07G024400','Bch09G052520','Bch08G025690','Bch01G001080','Bch08G025700','Bch09G022700','Bch02G037440','Bch06G049170','Bch03G052780','Bch01G011490')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # gSUS
    DotPlot(sce, group.by = group, features = c('Bch07G016400','Bch07G032990','Bch07G039810','Bch08G031040','Bch09G036340','Bch02G020530','Bch07G034310','Bch02G014380','Bch09G073230','Bch10G003630','Bch06G050090','Bch02G036340','Bch08G039940','Bch10G004100','Bch09G072990','Bch05G012500','Bch02G025520','Bch07G043080')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # hEP or hMCE or hPEN (shared marker)
    DotPlot(sce, group.by = group, features = c('Bch06G047330')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  }
  
  ## 3e. Cell cycle markers
  if(TRUE){
    # G0
    DotPlot(sce, group.by = group, features = c('Bch10G010480','Bch03G016810')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # G1
    DotPlot(sce, group.by = group, features = c('Bch06G046480','Bch06G046510','Bch06G046520','Bch06G046500','Bch06G046540','Bch09G008720')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # G1/S
    DotPlot(sce, group.by = group, features = c('Bch03G041140','Bch03G046750','Bch05G040980','Bch01G044550','Bch03G038030')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # S
    DotPlot(sce, group.by = group, features = c('Bch03G026510','Bch09G073880','Bch08G040370','Bch10G007910','Bch05G023080')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # G2
    DotPlot(sce, group.by = group, features = c('Bch03G021730','Bch09G064180','Bch04G006460')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # M
    DotPlot(sce, group.by = group, features = c('Bch09G030240','Bch07G044520','Bch07G017650','Bch06G022310','Bch03G065490','Bch06G035910','Bch09G011990')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  }
  
  ## 3f. FISH-validated markers
  if(TRUE){
    # FISH_Embryo
    DotPlot(sce, group.by = group, features = c('Bch01G037870','Bch05G031720','Bch02G013660','Bch02G007660','Bch03G004570','Bch02G003900','Bch06G034290','Bch06G002940','Bch05G020720','Bch03G067230','Bch01G001270','Bch08G025520','Bch10G030600','Bch10G025490','Bch05G011150','Bch04G029010','Bch09G040240','Bch03G025520','Bch05G005580','Bch06G002280','Bch08G005810','Bch0g005000','Bch04G005670','Bch07G026150','Bch07G013890','Bch04G036020','Bch07G025680','Bch03G025380','Bch08G037960','Bch09G070930','Bch03G027320','Bch09G069220','Bch08G007410','Bch03G046480','Bch02G040540','Bch09G017340','Bch07G032990','Bch07G039810','Bch07G016400','Bch08G031040','Bch07G041900','Bch02G024660','Bch06G012600','Bch09G066100','Bch07G028700')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_chalazal_cyst
    DotPlot(sce, group.by = group, features = c('Bch05G009200','Bch03G034530','Bch03G020620','Bch09G033990','Bch09G031150','Bch06G032040','Bch02G051350','Bch09G008670','Bch05G022740')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_chalazal_nodules
    DotPlot(sce, group.by = group, features = c('Bch08G000380')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_PEN, MCE, embryo
    DotPlot(sce, group.by = group, features = c('Bch09G028790')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_MCE, SC
    DotPlot(sce, group.by = group, features = c('Bch06G006810','Bch09G070540')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_ESR
    DotPlot(sce, group.by = group, features = c('Bch03G070800','Bch06G044300')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_oi2
    DotPlot(sce, group.by = group, features = c('Bch07G028300')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_ii1_oi1
    DotPlot(sce, group.by = group, features = c('Bch06G041480','Bch02G046020')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_ii1
    DotPlot(sce, group.by = group, features = c('Bch07G026360','Bch09G000340','Bch02G030100','Bch10G001340','Bch09G074580')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_ii1/ii2
    DotPlot(sce, group.by = group, features = c('Bch07G016590','Bch09G046750')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
    # FISH_oi1/oi2 (empty feature list – left as placeholder)
    DotPlot(sce, group.by = group, features = c('')) + 
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))
  }
}
