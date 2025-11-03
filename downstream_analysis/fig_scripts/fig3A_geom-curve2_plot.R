library(tidyverse)
library(readxl)
library(plyr)
library(vctrs)
library(ggplot2)
library(ggpubr)
library(grid)
library(ggsci)
source("geom-curve2.R")

show_col(pal_npg('nrc')(50))
show_col(pal_igv('default')(50))
show_col(pal_igv('alternating')(50))
setwd("scripts/")

metadata = read.table('../CK_batch3_10samples/CK3D/CK3D_metadata_res1.txt')
metadata = read.table('../CK_batch3_10samples/CK5D_subcluster/CK5D_metadata_res1_subcluster.txt')
metadata = read.table('../CK_batch3_10samples/CK7D/CK7D_metadata_res1.txt')
metadata = read.table('../CK_batch3_10samples/CK9D/CK9D_metadata_res1.txt')
metadata = read.table('../CK_batch3_10samples/CK11D_subcluster/CK11D_metadata_res1_subcluster.txt')
table(metadata$harmony_clusters_res1)
sum(table(metadata$harmony_clusters_res1))
a = table(metadata$harmony_clusters_res1)/sum(table(metadata$harmony_clusters_res1))
a[c(18,20,3,5,19,13,15,17,6,9,2,4,8,10,14,16,11,7,1,12)]
a[c(14,16,3,11,6,17,13,18,7,9,4,12,2,15,10,8,1,5)]
a[c(16,19,3,25,11,17,24,8,9,22,10,12,13,2,20,23,7,15,4,6,21,14,5,1,18)]
a[c(2,23,19,11,15,20,9,16,21,10,14,4,22,17,18,3,7,8,1,6,5,12,13)]
a[c(12,19,23,9,1,16,2,22,6,15,21,24,4,7,8,13,17,18,20,10,14,25,26,27,3,11,5)]

dat01 <- read_excel("../CK_batch3_10samples/CK_clusters.xlsx")
dat01$newX = dat01$newX - 0.5
dat02 <- read_excel("../CK_batch3_10samples/CK_clusters_cor.xlsx")
dat02$startx = dat02$startx - 0.5
dat02$endx = dat02$endx - 0.5

dat02 = dat02 %>% arrange(desc(size2))
top <- dat02 %>% group_by(startx,starty) %>% slice_head(n = 1) %>% ungroup()

ggplot() +
  geom_curve2(data = top, alpha = 0.5,
              aes(x = startx, xend=endx,
                  y = starty, yend=endy,
                  curvature = curvature,
                  size = size,
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
ggsave('../CK_batch3_10samples/CK_clusters_link.svg', plot = last_plot(), height = 8, width = 12)

# TT
if(TRUE){
  metadata = read.table('../TT_batch3_8samples/TT5D/TT5D_metadata_res1.txt')
  metadata = read.table('../TT_batch3_8samples/TT11D/TT11D_metadata_res1.txt')
  a = table(metadata$unintegrated_clusters_res1)/sum(table(metadata$unintegrated_clusters_res1))
  a[c(16,15,1,14,18,10,12,13,17,3,7,11,2,5,9,6,8,4)]
  a[c(10,13,15,3,1,11,12,14,16,17,2,4,5,6,8,9,7)]
  
  metadata = read.table('../TT_batch3_8samples/TT7D/TT7D_metadata_res1.txt')
  metadata = read.table('../TT_batch3_8samples/TT8D/TT8D_metadata_res1_rep1_rep3.txt')
  metadata = read.table('../TT_batch3_8samples/TT9D/TT9D_metadata_res1.txt')
  a = table(metadata$harmony_clusters_res1)/sum(table(metadata$harmony_clusters_res1))
  a[c(11,20,17,19,3,18,1,13,15,16,4,8,9,10,6,14,5,7,12,2)]
  a[c(20,14,16,12,13,2,23,19,1,15,17,18,21,22,3,7,9,8,5,10,4,11,6)]
  a[c(19,21,1,15,3,20,10,12,13,14,17,2,22,4,5,7,8,18,16,6,11,9)]
  
  dat01 <- read_excel("../TT_batch3_8samples/TT_cluster.xlsx")
  dat01$newX = dat01$newX - 0.5
  dat02 <- read_excel("../TT_batch3_8samples/TT_cluster_cor.xlsx")
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
  ggsave('../TT_batch3_8samples/TT_clusters_link.svg', plot = last_plot(), height = 8, width = 12)
  

}
