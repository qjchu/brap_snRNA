library(tidyverse)
library(readxl)
library(plyr)
library(vctrs)
library(ggplot2)
library(ggpubr)
library(grid)
library(ggsci)
library(sva)

show_col(pal_npg('nrc')(50))
show_col(pal_igv('default')(50))
show_col(pal_igv('alternating')(50))
setwd('D:/CHU THINKBOOK/A03 白菜修改/scripts/')
source("fig/geom-curve2.R")

# TT
if(TRUE){
  metadata = read.table('../revised_data/TT5D/TT5D_combined_metadata_res1_Louvain.txt')
  metadata = read.table('../revised_data/TT11D/TT11D_combined_metadata_res1.txt')
  a = table(metadata$unintegrated_clusters_res1)/sum(table(metadata$unintegrated_clusters_res1))
  # 根据excel里面的顺序排序
  a[c(18,24,15,8,10,14,17,23,5,11,19,20,2,6,13,4,3,16,1,7,9,12,21,22)] # TT5D
  a[c(17,1,11,9,13,5,2,3,4,6,7,8,10,12,14,15,16)] # TT11D
  
  metadata = read.table('../revised_data/TT7D/TT7D_combined_metadata_res1_5.txt')
  metadata = read.table('../revised_data/TT8D/TT8D_combined_metadata_res1_rep1_rep3.txt')
  metadata = read.table('../revised_data/TT9D/TT9D_combined_metadata_res1_Louvain.txt')
  a = table(metadata$harmony_clusters_res1)/sum(table(metadata$harmony_clusters_res1))
  # 根据excel里面的顺序排序
  a[c(11,25,22,16,17,19,23,9,4,14,2,20,6,1,10,21,24,3,5,7,8,12,13,15,18)] # TT7D
  a[c(15,19,22,24,17,12,2,9,21,7,6,13,4,10,5,23,1,3,8,11,14,16,18,20)] # TT8D
  a[c(20,28,24,12,13,15,21,5,22,19,4,10,18,7,8,9,26,1,2,3,6,11,14,16,17,23,25,27)] # TT9D
  
  # corr
  if(TRUE){
    exprSet21 <- read.table('../revised_data/TT5D/TT5D_combined_RNA_clusters_res1_Louvain_avg.txt', header = TRUE)
    exprSet22 <- read.table('../revised_data/TT7D/TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt', header = TRUE)
    exprSet23 <- read.table('../revised_data/TT8D/TT8D_combined_RNA_harmony_clusters_res1_avg_rep1_rep3.txt', header = TRUE)
    exprSet24 <- read.table('../revised_data/TT9D/TT9D_combined_RNA_harmony_clusters_res1_Louvain_avg.txt', header = TRUE)
    exprSet25 <- read.table('../revised_data/TT11D/TT11D_combined_RNA_clusters_res1_avg.txt', header = TRUE)
    exprSet2 = cbind(exprSet21, exprSet22[rownames(exprSet21),], exprSet23[rownames(exprSet21),], exprSet24[rownames(exprSet21),], exprSet25[rownames(exprSet21),])
    colnames(exprSet2) = c('TT5D SC (c1)', 'TT5D SC-endosperm (c2)', 'TT5D SC-oi (c3)', 'TT5D SC-ii (c4)', 'TT5D SC-dividing-cell (c5)', 'TT5D SC-endosperm (c6)', 'TT5D SC (c7)', 'TT5D CZSC (c8)', 'TT5D SC (c9)', 'TT5D CZSC (c10)', 'TT5D SC-dividing-cell (c11)', 'TT5D SC (c12)', 'TT5D SC-ii (c13)', 'TT5D CZSC (c14)', 'TT5D EP-dividing-cell (c15)', 'TT5D SC-SUS (c16)', 'TT5D CZSC (c17)', 'TT5D Endosperm-PCD (c18)', 'TT5D SC-dividing-cell (c19)', 'TT5D SC-dividing-cell (c20)', 'TT5D SC (c21)', 'TT5D SC (c22)', 'TT5D CZSC-phloem-xylem (c23)', 'TT5D Endosperm-active (c24)',
                           'TT7D SC-oi (c1)', 'TT7D SC-ii (c2)', 'TT7D SC (c3)', 'TT7D SC-endosperm (c4)', 'TT7D SC (c5)', 'TT7D SC-ii (c6)', 'TT7D SC (c7)', 'TT7D SC (c8)', 'TT7D SC-dividing-cell (c9)', 'TT7D SC-oi (c10)', 'TT7D Endosperm-PCD (c11)', 'TT7D SC (c12)', 'TT7D SC (c13)', 'TT7D SC-ii (c14)', 'TT7D SC (c15)', 'TT7D CZSC (c16)', 'TT7D CZSC (c17)', 'TT7D SC (c18)', 'TT7D CZSC (c19)', 'TT7D SC-ii (c20)', 'TT7D SC-oi (c21)', 'TT7D EP-dividing-cell (c22)', 'TT7D CZSC-phloem-xylem (c23)', 'TT7D SC-SUS (c24)', 'TT7D Endosperm-active (c25)',
                           'TT8D SC (c1)', 'TT8D CZSC (c2)', 'TT8D SC (c3)', 'TT8D SC-ii (c4)', 'TT8D SC-oi (c5)', 'TT8D SC-endosperm (c6)', 'TT8D SC-dividing-cell (c7)', 'TT8D SC (c8)', 'TT8D CZSC (c9)', 'TT8D SC-oi (c10)', 'TT8D SC (c11)', 'TT8D CZSC (c12)', 'TT8D SC-ii (c13)', 'TT8D SC (c14)', 'TT8D Endosperm-PCD (c15)', 'TT8D SC (c16)', 'TT8D EP-dividing-cell (c17)', 'TT8D SC (c18)', 'TT8D Endosperm-PCD (c19)', 'TT8D SC (c20)', 'TT8D CZSC-phloem-xylem (c21)', 'TT8D Endosperm-active (c22)', 'TT8D SC-SUS (c23)', 'TT8D Endosperm-active (c24)',
                           'TT9D SC (c1)', 'TT9D SC (c2)', 'TT9D SC (c3)', 'TT9D SC-endosperm (c4)', 'TT9D CZSC (c5)', 'TT9D SC (c6)', 'TT9D SC-ii (c7)', 'TT9D SC-oi (c8)', 'TT9D SC-oi (c9)', 'TT9D SC-ii (c10)', 'TT9D SC (c11)', 'TT9D CZSC (c12)', 'TT9D CZSC (c13)', 'TT9D SC (c14)', 'TT9D CZSC (c15)', 'TT9D SC (c16)', 'TT9D SC (c17)', 'TT9D SC-ii (c18)', 'TT9D SC-dividing-cell (c19)', 'TT9D Endosperm-PCD (c20)', 'TT9D CZSC (c21)', 'TT9D CZSC-phloem-xylem (c22)', 'TT9D SC (c23)', 'TT9D EP-dividing-cell (c24)', 'TT9D SC (c25)', 'TT9D SC-SUS (c26)', 'TT9D SC (c27)', 'TT9D Endosperm-active (c28)',
                           'TT11D CZSC (c1)', 'TT11D SC (c2)', 'TT11D SC (c3)', 'TT11D SC (c4)', 'TT11D SC-oi (c5)', 'TT11D SC (c6)', 'TT11D SC (c7)', 'TT11D SC (c8)', 'TT11D CZSC (c9)', 'TT11D SC (c10)', 'TT11D CZSC (c11)', 'TT11D SC (c12)', 'TT11D CZSC-phloem-xylem (c13)', 'TT11D SC (c14)', 'TT11D SC (c15)', 'TT11D SC (c16)', 'TT11D Endosperm (c17)')
    
    rownames(exprSet2)
    exprSet2 = exprSet2[rowSums(exprSet2) > 1,]
    exprSet2 = na.omit(exprSet2)
    
    max(exprSet2)
    min(exprSet2)
    
    dat_temp = exprSet2

    # remove batch -- combat
    if(TRUE){
      a = ComBat(scale(dat_temp), batch = c(rep('batch1',24), rep('batch2',25), rep('batch3',24), rep('batch4',28), rep('batch5',17)))
      
      # cor
      M3 <- cor(a, method = 'spearman')
      pheatmap::pheatmap(M3)
      ggsave('../revised_data/TT_cor.svg', plot = pheatmap::pheatmap(M3), height = 17, width = 18)
      
      ## TT5D TT7D corr
      zz = c()
      for (i in c(18,24,15,8,10,14,17,23,5,11,19,20,2,6,13,4,3,16,1,7,9,12,21,22)) {
        zz = c(zz, M3[i, 24 + c(11,25,22,16,17,19,23,9,4,14,2,20,6,1,10,21,24,3,5,7,8,12,13,15,18)])
      }
      write.table(zz, file = 'temp.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)

      ## TT7D TT8D corr
      zz = c()
      for (i in c(11,25,22,16,17,19,23,9,4,14,2,20,6,1,10,21,24,3,5,7,8,12,13,15,18)) {
        zz = c(zz, M3[24 + i, 24 + 25 + c(15,19,22,24,17,12,2,9,21,7,6,13,4,10,5,23,1,3,8,11,14,16,18,20)])
      }
      write.table(zz, file = 'temp.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
      
      ## TT8D TT9D corr
      zz = c()
      for (i in c(15,19,22,24,17,12,2,9,21,7,6,13,4,10,5,23,1,3,8,11,14,16,18,20)) {
        zz = c(zz, M3[24 + 25 + i, 24 + 25 + 24 + c(20,28,24,12,13,15,21,5,22,19,4,10,18,7,8,9,26,1,2,3,6,11,14,16,17,23,25,27)])
      }
      write.table(zz, file = 'temp.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
      
      ## TT9D TT11D corr
      zz = c()
      for (i in c(20,28,24,12,13,15,21,5,22,19,4,10,18,7,8,9,26,1,2,3,6,11,14,16,17,23,25,27)) {
        zz = c(zz, M3[24 + 25 + 24 + i, 24 + 25 + 24 + 28 + c(17,1,11,9,13,5,2,3,4,6,7,8,10,12,14,15,16)])
      }
      write.table(zz, file = 'temp.txt', quote = FALSE, col.names = FALSE, row.names = FALSE)
      
    }
  }
  
  dat01 <- read_excel("fig/fig3B_TT_cluster.xlsx")
  dat01$newX = dat01$newX - 0.5
  dat02 <- read_excel("fig/fig3B_TT_cluster_cor.xlsx")
  dat02$startx = dat02$startx - 0.5
  dat02$endx = dat02$endx - 0.5
  
  dat02 = dat02 %>% arrange(desc(size2))
  top <- dat02 %>% group_by(startx,starty) %>% slice_head(n = 1) %>% ungroup()
  
  ggplot() +
    geom_curve2(data = top, alpha = 0.5,
                aes(x = startx, xend=endx,
                    y = starty, yend=endy,
                    curvature = curvature,
                    size = size2,
                    color = group), 
                node.color = NA, ncp = 1, node.fill = NA) +
    scale_color_npg() +
    ggnewscale::new_scale_color() +
    ggnewscale::new_scale("size") +
    geom_point(data=dat01, aes(x = newX, y = newY, size = value3, color = group))+
    scale_size_continuous(range = c(0.5,15)) +
    scale_linewidth_continuous(range=c(0,0.1)) +
    scale_color_manual(values = c("Endosperm"="#e64b35", "EP"="#4dbbd5", "SC"="#00a087"))+
    annotate(geom = "label",x=dat01$newX, y = dat01$newY,label = dat01$celltype3, label.size = NA, fill="grey95", size = 3, alpha = 0.5) +
    theme_void() +
    xlim(0,5)
  ggsave('../revised_data/fig3B_TT_clusters_link.svg', plot = last_plot(), height = 8, width = 12)
  
}