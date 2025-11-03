library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)

getwd()
setwd("scripts")

# nFeature_RNA plot --> getDATA
if(TRUE){
  # CK
  rds = readRDS('../CK_batch3_10samples/10samples/CK_10samples_batch3_metadata_res1.rds')
  head(rds@meta.data)
  rds@meta.data$sample = gsub(pattern = '_.*', replacement = '', rds@meta.data$orig.ident)
  rds@meta.data$sample = factor(rds@meta.data$sample, levels = rev(c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')))
  vplot_data = VlnPlot(rds, features = c("nFeature_RNA"), ncol = 1, group.by ="sample", pt.size = 0)
  head(vplot_data[[1]][["data"]])
  write.table(x = vplot_data[[1]][["data"]], file = '../CK_batch3_10samples/CK_vlnplot_data_nFeature.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)

  vplot_data = VlnPlot(rds, features = c("nCount_RNA"), ncol = 1, group.by ="sample", pt.size = 0)
  write.table(x = vplot_data[[1]][["data"]], file = '../CK_batch3_10samples/CK_vlnplot_data_nCount.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)

  # TT
  rds = readRDS('../revised_data/TT_8samples/TT_8samples_combined_batch3_metadata_res1.rds')
  head(rds@meta.data)
  rds@meta.data$sample = gsub(pattern = '_.*', replacement = '', rds@meta.data$orig.ident)
  rds@meta.data$sample = factor(rds@meta.data$sample, levels = rev(c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')))
  vplot_data = VlnPlot(rds, features = c("nFeature_RNA"), ncol = 1, group.by ="sample", pt.size = 0)
  head(vplot_data[[1]][["data"]])
  write.table(x = vplot_data[[1]][["data"]], file = '../revised_data/TT_8samples/TT_combined_vlnplot_data_nFeature.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
  vplot_data = VlnPlot(rds, features = c("nCount_RNA"), ncol = 1, group.by ="sample", pt.size = 0)
  write.table(x = vplot_data[[1]][["data"]], file = '../revised_data/TT_8samples/TT_combined_vlnplot_data_nCount.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
  
}

# CK_vlnplot_data_nFeature.txt
vplot_data = read.table('../CK_batch3_10samples/10samples/CK_vlnplot_data_nFeature.txt')
vplot_data$ident = factor(vplot_data$ident, levels = rev(c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')))
quantile(vplot_data[vplot_data$ident == 'CK3D',1])
quantile(vplot_data[vplot_data$ident == 'CK5D',1])
quantile(vplot_data[vplot_data$ident == 'CK7D',1])
quantile(vplot_data[vplot_data$ident == 'CK9D',1])
quantile(vplot_data[vplot_data$ident == 'CK11D',1])
length(rownames(vplot_data[vplot_data$ident == 'CK3D',]))
length(rownames(vplot_data[vplot_data$ident == 'CK5D',]))
length(rownames(vplot_data[vplot_data$ident == 'CK7D',]))
length(rownames(vplot_data[vplot_data$ident == 'CK9D',]))
length(rownames(vplot_data[vplot_data$ident == 'CK11D',]))

ggplot(vplot_data, aes(factor(ident), nFeature_RNA, color = ident)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, cex=1.2) + geom_boxplot(width=0.1, cex=1.2) +
  scale_color_manual(values = rev(c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22'))) +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_blank(),
        axis.line = element_line(colour = 'grey60'),
        legend.position = 'none') 
ggsave('../CK_batch3_10samples/CK_vlnplot_nFeature_RNA2.svg', plot = last_plot(), height = 4, width = 1.6)

# CK_vlnplot_data_nCount.txt
vplot_data = read.table('../CK_batch3_10samples/10samples/CK_vlnplot_data_nCount.txt')
vplot_data$ident = factor(vplot_data$ident, levels = rev(c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')))
quantile(vplot_data[vplot_data$ident == 'CK3D',1])
quantile(vplot_data[vplot_data$ident == 'CK5D',1])
quantile(vplot_data[vplot_data$ident == 'CK7D',1])
quantile(vplot_data[vplot_data$ident == 'CK9D',1])
quantile(vplot_data[vplot_data$ident == 'CK11D',1])

ggplot(vplot_data, aes(factor(ident), nCount_RNA, color = ident)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, cex=1.2) + geom_boxplot(width=0.1, cex=1.2) +
  scale_color_manual(values = rev(c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22'))) +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_blank(),
        axis.line = element_line(colour = 'grey60'),
        legend.position = 'none') + 
  ylim(0,10000)
ggsave('../CK_batch3_10samples/CK_vlnplot_nCount_RNA2.svg', plot = last_plot(), height = 4, width = 1.6)


# TT nFeature_RNA
vplot_data = read.table('../revised_data/TT_8samples/TT_combined_vlnplot_data_nFeature.txt')
vplot_data$ident = factor(vplot_data$ident, levels = rev(c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')))
quantile(vplot_data[vplot_data$ident == 'TT5D',1])
quantile(vplot_data[vplot_data$ident == 'TT7D',1])
quantile(vplot_data[vplot_data$ident == 'TT8D',1])
quantile(vplot_data[vplot_data$ident == 'TT9D',1])
quantile(vplot_data[vplot_data$ident == 'TT11D',1])
length(rownames(vplot_data[vplot_data$ident == 'TT5D',]))
length(rownames(vplot_data[vplot_data$ident == 'TT7D',]))
length(rownames(vplot_data[vplot_data$ident == 'TT8D',]))
length(rownames(vplot_data[vplot_data$ident == 'TT9D',]))
length(rownames(vplot_data[vplot_data$ident == 'TT11D',]))

ggplot(vplot_data, aes(factor(ident), nFeature_RNA, color = ident)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, cex=1.2) + geom_boxplot(width=0.1, cex=1.2) +
  scale_color_manual(values = rev(c( '#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B'))) +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_blank(),
        axis.line = element_line(colour = 'grey60'),
        legend.position = 'none') 
ggsave('../revised_data/TT_8samples/TT_combined_vlnplot_nFeature.svg', plot = last_plot(), height = 4, width = 1.6)

# TT nCount
vplot_data = read.table('../revised_data/TT_8samples/TT_combined_vlnplot_data_nCount.txt')
vplot_data$ident = factor(vplot_data$ident, levels = rev(c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')))
quantile(vplot_data[vplot_data$ident == 'TT5D',1])
quantile(vplot_data[vplot_data$ident == 'TT7D',1])
quantile(vplot_data[vplot_data$ident == 'TT8D',1])
quantile(vplot_data[vplot_data$ident == 'TT9D',1])
quantile(vplot_data[vplot_data$ident == 'TT11D',1])
ggplot(vplot_data, aes(factor(ident), nCount_RNA, color = ident)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, cex=1.2) + geom_boxplot(width=0.1, cex=1.2) +
  scale_color_manual(values = rev(c( '#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B'))) +
  coord_flip() +
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90,vjust = 0.5,hjust = 0.5),
        axis.title = element_blank(),
        axis.line = element_line(colour = 'grey60'),
        legend.position = 'none') + 
  ylim(0,10000)
ggsave('../revised_data/TT_8samples/TT_combined_vlnplot_nCount.svg', plot = last_plot(), height = 4, width = 1.6)

