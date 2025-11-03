library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(ggsci)
library(ggrepel)
library(tidyr)
library(stringr)
library(sva)

sessionInfo()
setwd("scripts/")

m1 = read.table('../CK_batch3_10samples/CK3D/CK3D_res1_harmony_cluster_markers_anno.txt', header = TRUE, sep = "\t")
colnames(m1) = c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene','atha_GeneID','atha_Symbol','computational_description','full_name')
m1 <- m1 %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
m1_top <- m1 %>% group_by(cluster) %>% slice_head(n = 100) %>%ungroup()

m2 = read.table('../CK_batch3_10samples/CK5D_subcluster/CK5D_res1_subcluster_harmony_cluster_markers_anno.txt', header = TRUE, sep = "\t")
colnames(m2) = c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene','atha_GeneID','atha_Symbol','computational_description','full_name')
m2 <- m2 %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
m2_top <- m2 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()

m3 = read.table('../CK_batch3_10samples/CK7D/CK7D_res1_harmony_cluster_markers_anno.txt', header = TRUE, sep = "\t", quote = "")
colnames(m3) = c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene','atha_GeneID','atha_Symbol','computational_description','full_name')
m3 <- m3 %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
m3_top <- m3 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()

m4 = read.table('../CK_batch3_10samples/CK9D/CK9D_res1_harmony_cluster_markers_anno.txt', header = TRUE, sep = "\t")
colnames(m4) = c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene','atha_GeneID','atha_Symbol','computational_description','full_name')
m4 <- m4 %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
m4_top <- m4 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()

m5 = read.table('../CK_batch3_10samples/CK11D_subcluster/CK11D_res1_subcluster_harmony_cluster_markers_anno.txt', header = TRUE, sep = "\t", quote = "")
colnames(m5) = c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene','atha_GeneID','atha_Symbol','computational_description','full_name')
m5 <- m5 %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
m5_top <- m5 %>% group_by(cluster) %>% slice_head(n = 100) %>% ungroup()

# EP_Diving_Cell
markers1 = Reduce(intersect, list(m1_top[m1_top$cluster == 3,]$gene, 
                                 m2_top[m2_top$cluster == 11,]$gene, 
                                 m3_top[m3_top$cluster == 11,]$gene, m3_top[m3_top$cluster == 17,]$gene,
                                 m4_top[m4_top$cluster == 10,]$gene))
markers1 = c('Bch01G042330','Bch02G014800','Bch06G006190')

# EP
markers2 = Reduce(intersect, list(m3_top[m3_top$cluster == 25,]$gene,
                                  m4_top[m4_top$cluster == 21,]$gene,
                                  m5_top[m5_top$cluster == 4,]$gene, m5_top[m5_top$cluster == 7,]$gene, m5_top[m5_top$cluster == 8,]$gene, m5_top[m5_top$cluster == 13,]$gene, m5_top[m5_top$cluster == 15,]$gene, m5_top[m5_top$cluster == 21,]$gene, m5_top[m5_top$cluster == 24,]$gene))
markers2 = Reduce(intersect, list(m3_top[m3_top$cluster == 25,]$gene,
                                  m4_top[m4_top$cluster == 21,]$gene,
                                  m5_top[m5_top$cluster == 21,]$gene))

# CZE
markers3 = Reduce(intersect, list(m1_top[m1_top$cluster == 20,]$gene, 
                                  m2_top[m2_top$cluster == 14,]$gene, 
                                  m3_top[m3_top$cluster == 16,]$gene,
                                  m4_top[m4_top$cluster == 2,]$gene))

# MCE
markers4 = Reduce(intersect, list(m1_top[m1_top$cluster == 18,]$gene, 
                                  m2_top[m2_top$cluster == 16,]$gene, 
                                  m3_top[m3_top$cluster == 19,]$gene))

# PEN
markers5 = Reduce(intersect, list(m1_top[m1_top$cluster == 18,]$gene, 
                                  m2_top[m2_top$cluster == 3,]$gene, 
                                  m3_top[m3_top$cluster == 3,]$gene,
                                  m4_top[m4_top$cluster == 9,]$gene))

# CZSC
markers6 = Reduce(intersect, list(m1_top[m1_top$cluster == 5,]$gene, 
                                  m2_top[m2_top$cluster == 6,]$gene, 
                                  m3_top[m3_top$cluster == 8,]$gene,
                                  m4_top[m4_top$cluster == 4,]$gene))

# CZSC-phloem-xylem
markers7 = Reduce(intersect, list(m1_top[m1_top$cluster == 19,]$gene, 
                                  m2_top[m2_top$cluster == 17,]$gene, 
                                  m3_top[m3_top$cluster == 22,]$gene,
                                  m4_top[m4_top$cluster == 22,]$gene,
                                  m5_top[m5_top$cluster == 18,]$gene))

# SC-ii
markers8 = Reduce(intersect, list(m1_top[m1_top$cluster == 7,]$gene, 
                                  m2_top[m2_top$cluster == 8,]$gene, 
                                  m3_top[m3_top$cluster == 5,]$gene,
                                  m4_top[m4_top$cluster == 5,]$gene,
                                  m5_top[m5_top$cluster == 5,]$gene))

# SC-oi
markers9 = Reduce(intersect, list(m1_top[m1_top$cluster == 1,]$gene, 
                                  m2_top[m2_top$cluster == 1,]$gene, 
                                  m3_top[m3_top$cluster == 1,]$gene,
                                  m4_top[m4_top$cluster == 12,]$gene,
                                  m5_top[m5_top$cluster == 14,]$gene))

# SC-SUS
markers10 = Reduce(intersect, list(m1_top[m1_top$cluster == 16,]$gene, 
                                  m2_top[m2_top$cluster == 15,]$gene, 
                                  m3_top[m3_top$cluster == 21,]$gene))


markers = c(markers1[1:3], markers2[1:3], markers3[1:3], markers4[1:3], markers5[1:3],
            markers6[1:3], markers7[1:3], markers8[1:3], markers9[1:3], markers10[1:3])
markers = c("Bch01G042330", "Bch02G014800", "Bch06G006190", # EP-dividing
            "Bch01G045760", "Bch01G047930", "Bch07G013960", # EP
            "Bch06G008080", "Bch06G009710", "Bch07G011060", # CZE
            "Bch03G067580", "Bch04G033350", "Bch05G008080", # MCE
            "Bch03G072300", "Bch04G014210", "Bch09G041620", # PEN
            "Bch05G007030", "Bch01G000920", "Bch09G050570", # CZSC
            "Bch09G027690", "Bch02G031690", "Bch03G063320", # CZSC-phloem-xylem
            "Bch01G013080", "Bch09G017230", "Bch07G022460", # SC-ii
            "Bch01G041200", "Bch07G015970", "Bch01G011550", # SC-oi
            "Bch03G029400", "Bch05G001280", "Bch03G065950") # SUS

# read RDS CK
if(TRUE){
  # CK3D
  sce = readRDS("../RDS/CK3D_merge_res1.rds")
  dotplot_data_CK3D <- DotPlot(sce, features = markers)
  dotplot_data_CK3D[["data"]]$id = paste('CK3D_', dotplot_data_CK3D[["data"]]$id, sep = '')
  dotplot_data_CK3D[["data"]]$time = c('CK3D')
  head(dotplot_data_CK3D[["data"]])
  # write.table(x = dotplot_data_CK3D[["data"]], file = '../CK_batch3_10samples/CK3D/CK3D_markers_dotplot_data.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # CK5D
  sce = readRDS("../RDS/CK5D_merge_res1_subcluster.rds")
  dotplot_data_CK5D <- DotPlot(sce, features = markers)
  dotplot_data_CK5D[["data"]]$id = paste('CK5D_', dotplot_data_CK5D[["data"]]$id, sep = '')
  dotplot_data_CK5D[["data"]]$time = c('CK5D')
  head(dotplot_data_CK5D[["data"]])
  # write.table(x = dotplot_data_CK5D[["data"]], file = '../CK_batch3_10samples/CK5D_subcluster/CK5D_markers_dotplot_data.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # CK7D
  sce = readRDS("../RDS/CK7D_merge_res1.rds")
  dotplot_data_CK7D <- DotPlot(sce, features = markers)
  dotplot_data_CK7D[["data"]]$id = paste('CK7D_', dotplot_data_CK7D[["data"]]$id, sep = '')
  dotplot_data_CK7D[["data"]]$time = c('CK7D')
  head(dotplot_data_CK7D[["data"]])
  # write.table(x = dotplot_data_CK7D[["data"]], file = '../CK_batch3_10samples/CK7D/CK7D_markers_dotplot_data.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # CK9D
  sce = readRDS("../RDS/CK9D_merge_res1.rds")
  dotplot_data_CK9D <- DotPlot(sce, features = markers)
  dotplot_data_CK9D[["data"]]$id = paste('CK9D_', dotplot_data_CK9D[["data"]]$id, sep = '')
  dotplot_data_CK9D[["data"]]$time = c('CK9D')
  head(dotplot_data_CK9D[["data"]])
  # write.table(x = dotplot_data_CK9D[["data"]], file = '../CK_batch3_10samples/CK9D/CK9D_markers_dotplot_data.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  # CK11D
  sce = readRDS("../RDS/CK11D_merge_res1_subcluster.rds")
  dotplot_data_CK11D <- DotPlot(sce, features = markers)
  dotplot_data_CK11D[["data"]]$id = paste('CK11D_', dotplot_data_CK11D[["data"]]$id, sep = '')
  dotplot_data_CK11D[["data"]]$time = c('CK11D')
  head(dotplot_data_CK11D[["data"]])
  # write.table(x = dotplot_data_CK11D[["data"]], file = '../CK_batch3_10samples/CK11D_subcluster/CK11D_markers_dotplot_data.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  dotplot_data = rbind(dotplot_data_CK3D[["data"]], dotplot_data_CK5D[["data"]], dotplot_data_CK7D[["data"]], dotplot_data_CK9D[["data"]], dotplot_data_CK11D[["data"]])
  
  table(dotplot_data$id)
  dotplot_data$id = factor(dotplot_data$id, levels = c('CK3D_3','CK3D_20','CK3D_18','CK3D_5','CK3D_19','CK3D_7','CK3D_11','CK3D_1','CK3D_12','CK3D_16',# 'CK3D_2','CK3D_4','CK3D_6','CK3D_8','CK3D_9','CK3D_10','CK3D_13','CK3D_14','CK3D_15','CK3D_17',
                                                       'CK5D_11','CK5D_14','CK5D_16','CK5D_3','CK5D_6', 'CK5D_17','CK5D_8','CK5D_10', 'CK5D_1','CK5D_5', 'CK5D_15',#'CK5D_2','CK5D_4','CK5D_7','CK5D_9','CK5D_12','CK5D_13','CK5D_18',
                                                       'CK7D_11','CK7D_17','CK7D_25','CK7D_16','CK7D_19','CK7D_3','CK7D_8','CK7D_9','CK7D_24','CK7D_22','CK7D_5','CK7D_14','CK7D_1','CK7D_18','CK7D_21',# 'CK7D_2','CK7D_4','CK7D_6','CK7D_7','CK7D_10','CK7D_12','CK7D_13','CK7D_15','CK7D_20','CK7D_23',
                                                       'CK9D_16','CK9D_21','CK9D_2','CK9D_23','CK9D_1','CK9D_9','CK9D_11','CK9D_15','CK9D_19','CK9D_20','CK9D_4','CK9D_14','CK9D_22','CK9D_5','CK9D_7','CK9D_12','CK9D_13',# 'CK9D_3','CK9D_6','CK9D_8','CK9D_10','CK9D_17','CK9D_18',
                                                       'CK11D_13','CK11D_4','CK11D_7','CK11D_8','CK11D_15', 'CK11D_21', 'CK11D_24','CK11D_9','CK11D_12','CK11D_19','CK11D_23','CK11D_1','CK11D_2','CK11D_6','CK11D_16','CK11D_22', 'CK11D_17','CK11D_20','CK11D_25','CK11D_18','CK11D_5','CK11D_26','CK11D_14'))#'CK11D_3','CK11D_10','CK11D_11','CK11D_27'))
  table(is.na(dotplot_data$id))
  dotplot_data = na.omit(dotplot_data)
  dotplot_data$features.plot = factor(dotplot_data$features.plot, levels = markers)
  dotplot_data$time = factor(dotplot_data$time, levels = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D'))
  
  dotplot_data[dotplot_data$features.plot %in% markers[1:3],'celltype'] = 'EP-dividing-cell'
  dotplot_data[dotplot_data$features.plot %in% markers[4:6],'celltype'] = 'EP'
  dotplot_data[dotplot_data$features.plot %in% markers[7:9],'celltype'] = 'CZE'
  dotplot_data[dotplot_data$features.plot %in% markers[10:12],'celltype'] = 'MCE'
  dotplot_data[dotplot_data$features.plot %in% markers[13:15],'celltype'] = 'PEN'
  dotplot_data[dotplot_data$features.plot %in% markers[16:18],'celltype'] = 'CZSC'
  dotplot_data[dotplot_data$features.plot %in% markers[19:21],'celltype'] = 'CZSC-phloem-xylem'
  dotplot_data[dotplot_data$features.plot %in% markers[22:24],'celltype'] = 'SC-ii'
  dotplot_data[dotplot_data$features.plot %in% markers[25:27],'celltype'] = 'SC-oi'
  dotplot_data[dotplot_data$features.plot %in% markers[28:30],'celltype'] = 'SC-SUS'
  dotplot_data$celltype = factor(dotplot_data$celltype, levels = rev(c('EP-dividing-cell', 'EP','CZE','MCE','PEN','CZSC','CZSC-phloem-xylem','SC-ii','SC-oi','SC-SUS')))
  
  ggplot(dotplot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
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
                          midpoint = mean(dotplot_data$avg.exp.scaled), 
                          limit = range(dotplot_data$avg.exp.scaled), 
                          space = "Lab", na.value = "grey50", guide = "colourbar", 
                          aesthetics = "color") + 
    facet_grid(dotplot_data$celltype~dotplot_data$time, scales = 'free', space = 'free') +
    coord_flip() +
    scale_y_discrete(labels = c('CK3D_1' = "SC-oi (c1)",  'CK3D_2' = "SC-dividing-cell (c2)",  'CK3D_3' = "EP-dividing-cell (c3)",  'CK3D_4' = "SC-dividing-cell (c4)",  'CK3D_5' = "CZSC (c5)", 
                                'CK3D_6' = "SC (c6)",  'CK3D_7' = "SC-ii (c7)",  'CK3D_8' = "SC-dividing-cell (c8)", 'CK3D_9' =  "SC (c9)",  'CK3D_10' = "SC_Endosperm (c10)",
                                'CK3D_11' = "SC-ii (c11)", 'CK3D_12' = "SC-oi (c12)", 'CK3D_13' = "SC (c13)", 'CK3D_14' = "SC_Endosperm (c14)", 'CK3D_15' = "SC (c15)",
                                'CK3D_16' = "SC-SUS (c16)", 'CK3D_17' = "SC (c17)", 'CK3D_18' = "PEN/MCE (c18)", 'CK3D_19' = "CZSC-phloem-xylem (c19)", 'CK3D_20' = "CZE (c20)",
                                'CK5D_1' = "SC-oi (c1)",  'CK5D_2' = "SC_Endosperm (c2)",  'CK5D_3' = "PEN (c3)", 
                                'CK5D_4' = "SC-dividing-cell (c4)",'CK5D_5' = "SC-oi (c5)",'CK5D_6' = "CZSC (c6)", 
                                'CK5D_7' = "SC (c7)", 'CK5D_8' = "SC-ii (c8)",'CK5D_9' =  "SC (c9)", 
                                'CK5D_10' = "SC-ii (c10)", 'CK5D_11' = "EP-dividing-cell (c11)", 'CK5D_12' = "SC_Endosperm (c12)",
                                'CK5D_13' = "SC (c13)", 'CK5D_14' = "CZE (c14)", 'CK5D_15' = "SC-SUS (c15)", 'CK5D_16' = "MCE (c16)", 'CK5D_17' = "CZSC-phloem-xylem (c17)", 'CK5D_18' = "SC (c18)",
                                'CK7D_1' = "SC-oi (c1)", 'CK7D_2' = "SC (c2)", 'CK7D_3' = "PEN (c3)", 'CK7D_4' = "SC_Endosperm (c4)", 'CK7D_5' = "SC-ii (c5)", 
                                'CK7D_6' = "SC_Endosperm (c6)", 'CK7D_7' = "SC-dividing-cell (c7)",
                                'CK7D_8' = "CZSC (c8)", 'CK7D_9' = "CZSC (c9)", 'CK7D_10' = "SC (c10)", 'CK7D_11' = "EP-dividing-cell (c11)", 'CK7D_12' = "SC (c12)",'CK7D_13' = "SC (c13)", 'CK7D_14' = "SC-ii (c14)",
                                'CK7D_15' = "SC_Endosperm (c15)", 'CK7D_16' = "CZE (c16)", 'CK7D_17' = "EP-dividing-cell (c17)", 'CK7D_18' = "SC-oi (c18)", 'CK7D_19' = "MCE (c19)", 'CK7D_20' = "SC (c20)", 'CK7D_21' = "SC-SUS (c21)",
                                'CK7D_22' = "CZSC-phloem-xylem (c22)", 'CK7D_23' = "SC (c23)", 'CK7D_24' = "CZSC (c24)", 'CK7D_25' = "EP (c25)",
                                'CK9D_1' = "PEN (c1)", 'CK9D_2' = "CZE (c2)", 'CK9D_3' = "SC (c3)", 'CK9D_4' = "CZSC (c4)", 'CK9D_5' = "SC-ii (c5)", 'CK9D_6' = "SC_Endosperm (c6)", 'CK9D_7' = "SC-oi (c7)", 'CK9D_8' = "SC (c8)",
                                'CK9D_9' = "PEN (c9)", 'CK9D_10' = "EP-dividing-cell (c10)", 'CK9D_11' = "PEN (c11)", 'CK9D_12' = "SC-oi (c12)", 'CK9D_13' = "SC-oi (c13)", 'CK9D_14' = "CZSC (c14)", 'CK9D_15' = "PEN (c15)", 'CK9D_16' = "PEN_around_EP (c16)",
                                'CK9D_17' = "SC (c17)", 'CK9D_18' = "SC (c18)", 'CK9D_19' = "MCE_PEN (c19)", 'CK9D_20' = "PEN (c20)", 'CK9D_21' = "EP (c21)", 'CK9D_22' = "CZSC-phloem-xylem (c22)", 'CK9D_23' = "CZE (c23)",
                                'CK11D_1' = "PEN (c1)",  'CK11D_2' = "PEN (c2)", 'CK11D_3' =  "SC (c3)", 
                                'CK11D_4' = "EP (c4)",  'CK11D_5' = "SC-ii (c5)", 'CK11D_6' = "PEN (c6)",  'CK11D_7' = "EP (c7)", 'CK11D_8' =  "EP (c8)", 
                                'CK11D_9' = "CZE (c9)", 'CK11D_10' =  "SC (c10)",
                                'CK11D_11' = "SC_PEN (c11)", 'CK11D_12' = "CZE (c12)",'CK11D_13' =  "EP-dividing-cell (c13)", 
                                'CK11D_14' = "SC-oi (c14)", 'CK11D_15' = "EP (c15)",
                                'CK11D_16' = "PEN (c16)", 'CK11D_17' = "CZSC (c17)", 'CK11D_18' = "CZSC-phloem-xylem (c18)", 
                                'CK11D_19' = "CZE (c19)", 'CK11D_20' = "CZSC (c20)",
                                'CK11D_21' = "EP (c21)", 'CK11D_22' = "PEN (c22)",'CK11D_23' =  "CZE (c23)", 
                                'CK11D_24' = "EP (c24)", 'CK11D_25' = "CZSC (c25)",
                                'CK11D_26' = "SC-ii (c26)", 'CK11D_27' = "SC (c27)"))
  ggsave('../revised_data/CK_markers_dotplot_top3.svg', plot = last_plot(), height = 7, width = 16)

}

# read RDS TT
if(TRUE){
  # TT5D
  sce = readRDS("../revised_data/TT5D/TT5D_combined_merge_res1_Louvain.rds")
  dotplot_data_CK3D <- DotPlot(sce, features = markers)
  dotplot_data_CK3D[["data"]]$id = paste('TT5D_', dotplot_data_CK3D[["data"]]$id, sep = '')
  dotplot_data_CK3D[["data"]]$time = c('TT5D')
  head(dotplot_data_CK3D[["data"]])

  # TT7D
  sce = readRDS("../revised_data/TT7D/TT7D_combined_merge_res1_5.rds")
  dotplot_data_CK5D <- DotPlot(sce, features = markers)
  dotplot_data_CK5D[["data"]]$id = paste('TT7D_', dotplot_data_CK5D[["data"]]$id, sep = '')
  dotplot_data_CK5D[["data"]]$time = c('TT7D')
  head(dotplot_data_CK5D[["data"]])

  # TT8D
  sce = readRDS("../revised_data/TT8D/TT8D_combined_merge_res1_rep1_rep3.rds")
  dotplot_data_CK7D <- DotPlot(sce, features = markers)
  dotplot_data_CK7D[["data"]]$id = paste('TT8D_', dotplot_data_CK7D[["data"]]$id, sep = '')
  dotplot_data_CK7D[["data"]]$time = c('TT8D')
  head(dotplot_data_CK7D[["data"]])

  # TT9D
  sce = readRDS("../revised_data/TT9D/TT9D_combined_merge_res1_Louvain.rds")
  dotplot_data_CK9D <- DotPlot(sce, features = markers)
  dotplot_data_CK9D[["data"]]$id = paste('TT9D_', dotplot_data_CK9D[["data"]]$id, sep = '')
  dotplot_data_CK9D[["data"]]$time = c('TT9D')
  head(dotplot_data_CK9D[["data"]])

  # TT11D
  sce = readRDS("../revised_data/TT11D/TT11D_combined_merge_res1.rds")
  dotplot_data_CK11D <- DotPlot(sce, features = markers)
  dotplot_data_CK11D[["data"]]$id = paste('TT11D_', dotplot_data_CK11D[["data"]]$id, sep = '')
  dotplot_data_CK11D[["data"]]$time = c('TT11D')
  head(dotplot_data_CK11D[["data"]])

  dotplot_data = rbind(dotplot_data_CK3D[["data"]], dotplot_data_CK5D[["data"]], dotplot_data_CK7D[["data"]], dotplot_data_CK9D[["data"]], dotplot_data_CK11D[["data"]])
  
  table(dotplot_data$id)
  dotplot_data$id = factor(dotplot_data$id, levels = c('TT5D_14','TT5D_23','TT5D_17','TT5D_7','TT5D_9','TT5D_13','TT5D_16','TT5D_22','TT5D_3','TT5D_12','TT5D_2','TT5D_15',
                                                       'TT7D_22','TT7D_25','TT7D_11','TT7D_16','TT7D_17','TT7D_19','TT7D_23','TT7D_2','TT7D_6','TT7D_14','TT7D_20','TT7D_','TT7D_1','TT7D_10','TT7D_21','TT7D_24',
                                                       'TT8D_17','TT8D_22','TT8D_24','TT8D_15','TT8D_19','TT8D_2','TT8D_9','TT8D_12','TT8D_21','TT8D_4','TT8D_13','TT8D_5','TT8D_10','TT8D_23',
                                                       'TT9D_23','TT9D_27','TT9D_19','TT9D_4','TT9D_11','TT9D_12','TT9D_20','TT9D_21','TT9D_6','TT9D_9','TT9D_17','TT9D_7','TT9D_8','TT9D_25',
                                                       'TT11D_17','TT11D_1','TT11D_9','TT11D_11','TT11D_13','TT11D_5'))
  table(is.na(dotplot_data$id))
  dotplot_data = na.omit(dotplot_data)
  dotplot_data$features.plot = factor(dotplot_data$features.plot, levels = markers)
  dotplot_data$time = factor(dotplot_data$time, levels = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D'))
  
  dotplot_data[dotplot_data$features.plot %in% markers[1:3],'celltype'] = 'EP-dividing-cell'
  dotplot_data[dotplot_data$features.plot %in% markers[4:6],'celltype'] = 'EP'
  dotplot_data[dotplot_data$features.plot %in% markers[7:9],'celltype'] = 'CZE'
  dotplot_data[dotplot_data$features.plot %in% markers[10:12],'celltype'] = 'MCE'
  dotplot_data[dotplot_data$features.plot %in% markers[13:15],'celltype'] = 'PEN'
  dotplot_data[dotplot_data$features.plot %in% markers[16:18],'celltype'] = 'CZSC'
  dotplot_data[dotplot_data$features.plot %in% markers[19:21],'celltype'] = 'CZSC-phloem-xylem'
  dotplot_data[dotplot_data$features.plot %in% markers[22:24],'celltype'] = 'SC-ii'
  dotplot_data[dotplot_data$features.plot %in% markers[25:27],'celltype'] = 'SC-oi'
  dotplot_data[dotplot_data$features.plot %in% markers[28:30],'celltype'] = 'SC-SUS'
  dotplot_data$celltype = factor(dotplot_data$celltype, levels = rev(c('EP-dividing-cell', 'EP','CZE','MCE','PEN','CZSC','CZSC-phloem-xylem','SC-ii','SC-oi','SC-SUS')))
  
  ggplot(dotplot_data, aes(x = features.plot, y = id, size = pct.exp, color = avg.exp.scaled)) +
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
                          midpoint = mean(dotplot_data$avg.exp.scaled), 
                          limit = range(dotplot_data$avg.exp.scaled), 
                          space = "Lab", na.value = "grey50", guide = "colourbar", 
                          aesthetics = "color") + 
    facet_grid(dotplot_data$celltype~dotplot_data$time, scales = 'free', space = 'free') +
    coord_flip() +
    scale_y_discrete(labels = c('TT5D_0' = 'SC (c1)', 'TT5D_1' = 'SC-endosperm (c2)', 'TT5D_2' = 'SC-oi (c3)', 'TT5D_3' = 'SC-ii (c4)', 'TT5D_4' = 'SC-dividing-cell (c5)', 'TT5D_5' = 'SC-endosperm (c6)', 'TT5D_6' = 'SC (c7)', 'TT5D_7' = 'CZSC (c8)', 'TT5D_8' = 'SC (c9)', 'TT5D_9' = 'CZSC (c10)', 'TT5D_10' = 'SC-dividing-cell (c11)', 'TT5D_11' = 'SC (c12)', 'TT5D_12' = 'SC-ii (c13)', 'TT5D_13' = 'CZSC (c14)', 'TT5D_14' = 'EP-dividing-cell (c15)', 'TT5D_15' = 'SC-SUS (c16)', 'TT5D_16' = 'CZSC (c17)', 'TT5D_17' = 'Endosperm-PCD (c18)', 'TT5D_18' = 'SC-dividing-cell (c19)', 'TT5D_19' = 'SC-dividing-cell (c20)', 'TT5D_20' = 'SC (c21)', 'TT5D_21' = 'SC (c22)', 'TT5D_22' = 'CZSC-phloem-xylem (c23)', 'TT5D_23' = 'Endosperm-active (c24)',
                                'TT7D_1' = 'SC-oi (c1)', 'TT7D_2' = 'SC-ii (c2)', 'TT7D_3' = 'SC (c3)', 'TT7D_4' = 'SC-endosperm (c4)', 'TT7D_5' = 'SC (c5)', 'TT7D_6' = 'SC-ii (c6)', 'TT7D_7' = 'SC (c7)', 'TT7D_8' = 'SC (c8)', 'TT7D_9' = 'SC-dividing-cell (c9)', 'TT7D_10' = 'SC-oi (c10)', 'TT7D_11' = 'Endosperm-PCD (c11)', 'TT7D_12' = 'SC (c12)', 'TT7D_13' = 'SC (c13)', 'TT7D_14' = 'SC-ii (c14)', 'TT7D_15' = 'SC (c15)', 'TT7D_16' = 'CZSC (c16)', 'TT7D_17' = 'CZSC (c17)', 'TT7D_18' = 'SC (c18)', 'TT7D_19' = 'CZSC (c19)', 'TT7D_20' = 'SC-ii (c20)', 'TT7D_21' = 'SC-oi (c21)', 'TT7D_22' = 'EP-dividing-cell (c22)', 'TT7D_23' = 'CZSC-phloem-xylem (c23)', 'TT7D_24' = 'SC-SUS (c24)', 'TT7D_25' = 'Endosperm-active (c25)',
                                'TT8D_1' = 'SC (c1)', 'TT8D_2' = 'CZSC (c2)', 'TT8D_3' = 'SC (c3)', 'TT8D_4' = 'SC-ii (c4)', 'TT8D_5' = 'SC-oi (c5)', 'TT8D_6' = 'SC-endosperm (c6)', 'TT8D_7' = 'SC-dividing-cell (c7)', 'TT8D_8' = 'SC (c8)', 'TT8D_9' = 'CZSC (c9)', 'TT8D_10' = 'SC-oi (c10)', 'TT8D_11' = 'SC (c11)', 'TT8D_12' = 'CZSC (c12)', 'TT8D_13' = 'SC-ii (c13)', 'TT8D_14' = 'SC (c14)', 'TT8D_15' = 'Endosperm-PCD (c15)', 'TT8D_16' = 'SC (c16)', 'TT8D_17' = 'EP-dividing-cell (c17)', 'TT8D_18' = 'SC (c18)', 'TT8D_19' = 'Endosperm-PCD (c19)', 'TT8D_20' = 'SC (c20)', 'TT8D_21' = 'CZSC-phloem-xylem (c21)', 'TT8D_22' = 'Endosperm-active (c22)', 'TT8D_23' = 'SC-SUS (c23)', 'TT8D_24' = 'Endosperm-active (c24)',
                                'TT9D_0' = 'SC (c1)', 'TT9D_1' = 'SC (c2)', 'TT9D_2' = 'SC (c3)', 'TT9D_3' = 'SC-endosperm (c4)', 'TT9D_4' = 'CZSC (c5)', 'TT9D_5' = 'SC (c6)', 'TT9D_6' = 'SC-ii (c7)', 'TT9D_7' = 'SC-oi (c8)', 'TT9D_8' = 'SC-oi (c9)', 'TT9D_9' = 'SC-ii (c10)', 'TT9D_10' = 'SC (c11)', 'TT9D_11' = 'CZSC (c12)', 'TT9D_12' = 'CZSC (c13)', 'TT9D_13' = 'SC (c14)', 'TT9D_14' = 'CZSC (c15)', 'TT9D_15' = 'SC (c16)', 'TT9D_16' = 'SC (c17)', 'TT9D_17' = 'SC-ii (c18)', 'TT9D_18' = 'SC-dividing-cell (c19)', 'TT9D_19' = 'Endosperm-PCD (c20)', 'TT9D_20' = 'CZSC (c21)', 'TT9D_21' = 'CZSC-phloem-xylem (c22)', 'TT9D_22' = 'SC (c23)', 'TT9D_23' = 'EP-dividing-cell (c24)', 'TT9D_24' = 'SC (c25)', 'TT9D_25' = 'SC-SUS (c26)', 'TT9D_26' = 'SC (c27)', 'TT9D_27' = 'Endosperm-active (c28)',
                                'TT11D_1' = 'CZSC (c1)', 'TT11D_2' = 'SC (c2)', 'TT11D_3' = 'SC (c3)', 'TT11D_4' = 'SC (c4)', 'TT11D_5' = 'SC-oi (c5)', 'TT11D_6' = 'SC (c6)', 'TT11D_7' = 'SC (c7)', 'TT11D_8' = 'SC (c8)', 'TT11D_9' = 'CZSC (c9)', 'TT11D_10' = 'SC (c10)', 'TT11D_11' = 'CZSC (c11)', 'TT11D_12' = 'SC (c12)', 'TT11D_13' = 'CZSC-phloem-xylem (c13)', 'TT11D_14' = 'SC (c14)', 'TT11D_15' = 'SC (c15)', 'TT11D_16' = 'SC (c16)', 'TT11D_17' = 'Endosperm (c17)'))
  ggsave('../revised_data/TT_8samples/TT_markers_dotplot_top3.svg', plot = last_plot(), height = 7, width = 16)

}

