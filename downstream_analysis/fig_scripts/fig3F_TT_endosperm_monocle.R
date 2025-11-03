library(monocle)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(ggsci)
library(ggrepel)
library(scales)
library(tidyr)
library(ggpubr)
library(clusterProfiler)
library(org.Brapaoleracea.eg.db)

show_col(pal_igv('default')(50))
setwd('scripts/')

cds = readRDS('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D/TT_combined_Endo_5D7D8D9D11D_monocle.rds')

plot_cell_trajectory(cds, color_by = 'Pseudotime', size=1, show_backbone=TRUE, alpha = 0.5) + 
  scale_color_bs5('red') +
  theme(legend.position = 'right')
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_time.png', plot = last_plot(), height = 4, width = 7, bg = "transparent")
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_time.svg', plot = last_plot(), height = 4, width = 7, bg = "transparent")

plot_cell_trajectory(cds, color_by = 'celltype', size=1, show_backbone=TRUE, alpha = 0.5) + 
  theme(legend.position = 'right') 
plot_cell_trajectory(cds, color_by = 'celltype', size=1, show_backbone=TRUE, alpha = 0.5) + 
  theme(legend.position = 'right') + facet_wrap("~celltype", nrow = 1)
plot_cell_trajectory(cds, color_by = 'State', size=1, show_backbone=TRUE, alpha = 0.5) +
  theme(legend.position = 'right')
cds = orderCells(cds, root_state = 2)


cds@phenoData@data[["celltype2"]] = cds@phenoData@data[["celltype"]]
cds@phenoData@data[["celltype2"]] = gsub('TT11D_Endo_c17', 'TT11D_Endosperm', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT5D_Endo_c18', 'TT5D_Endosperm_PCD', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT5D_Endo_c24', 'TT5D_Endosperm_active', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT7D_Endo_c11', 'TT7D_Endosperm_PCD', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT7D_Endo_c25', 'TT7D_Endosperm_active', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT8D_Endo_c15', 'TT8D_Endosperm_PCD', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT8D_Endo_c19', 'TT8D_Endosperm_PCD', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT8D_Endo_c22', 'TT8D_Endosperm_active', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT8D_Endo_c24', 'TT8D_Endosperm_active', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT9D_Endo_c20', 'TT9D_Endosperm_PCD', cds@phenoData@data[["celltype2"]])
cds@phenoData@data[["celltype2"]] = gsub('TT9D_Endo_c28', 'TT9D_Endosperm_active', cds@phenoData@data[["celltype2"]])
table(cds@phenoData@data[["celltype2"]])
plot_cell_trajectory(cds, color_by = 'celltype2', size=1, show_backbone=TRUE, alpha = 0.5) + 
  theme(legend.position = 'right') 
cds@phenoData@data[["celltype3"]] = cds@phenoData@data[["celltype2"]]
cds@phenoData@data[["celltype3"]] = gsub('TT11D_Endosperm', 'Endosperm_PCD', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT5D_Endosperm_PCD', 'Endosperm_PCD', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT7D_Endosperm_PCD', 'Endosperm_PCD', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT8D_Endosperm_PCD', 'Endosperm_PCD', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT9D_Endosperm_PCD', 'Endosperm_PCD', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT5D_Endosperm_active', 'Endosperm_active', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT7D_Endosperm_active', 'Endosperm_active', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT8D_Endosperm_active', 'Endosperm_active', cds@phenoData@data[["celltype3"]])
cds@phenoData@data[["celltype3"]] = gsub('TT9D_Endosperm_active', 'Endosperm_active', cds@phenoData@data[["celltype3"]])
plot_cell_trajectory(cds, color_by = 'celltype3', size=1, show_backbone=TRUE, alpha = 0.5) + 
  theme(legend.position = 'right') + scale_color_igv()
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D/TT_combined_Endo_5D7D8D9D11D_monocle_type.svg', plot = last_plot(), height = 4, width = 7, bg = "transparent")


cds@phenoData@data$time = gsub('_.*', '', cds@phenoData@data$celltype)
cds@phenoData@data$time = factor(cds@phenoData@data$time, levels = c('TT5D', 'TT7D', 'TT8D', 'TT9D','TT11D'))
plot_cell_trajectory(cds, color_by = 'time', size=1, show_backbone=TRUE) + 
  scale_color_manual(values = c('#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
  theme(legend.position = 'right')
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle.svg', plot = last_plot(), height = 4, width = 7, bg = "transparent")

plot_cell_trajectory(cds, color_by = 'time', size=1, show_backbone=TRUE) + 
  scale_color_manual(values = c('#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
  theme(legend.position = 'right') + facet_wrap("~time", nrow = 1)
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_facet.png', plot = last_plot(), height = 4, width = 15, bg = "transparent")
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_facet.svg', plot = last_plot(), height = 4, width = 15, bg = "transparent")
plot_cell_trajectory(cds, color_by = 'celltype', size=1, show_backbone=TRUE) + 
  theme(legend.position = 'right') + facet_wrap("~time", nrow = 1)
plot_cell_trajectory(cds, color_by = 'Pseudotime', size=1, show_backbone=TRUE) + 
  scale_color_bs5('red') +
  theme(legend.position = 'right') + facet_wrap("~time", nrow = 1)

ggplot(pData(cds), aes(Pseudotime, colour = time, fill=time)) +
  scale_color_manual(values = c('#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
  scale_fill_manual(values = c('#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B')) +
  geom_density(bw=0.5,size=1,alpha = 0.5) + theme_classic2()
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_density.svg', plot = last_plot(), height = 4, width = 7, bg = "transparent")



deg = read.table('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_DEG.txt')

keygenes = c('Bch02G005190','Bch09G046750','Bo2g011110','Bch06G036690','Bch06G019890','Bch06G029280') 
# State 3
keygenes = c('Bch06G010520','Bch07G028630','Bch09G069080','Bch01G050010','Bch05G049660','Bch02G027190','Bch07G045250') 
plot_genes_in_pseudotime(cds[keygenes,], color_by = "celltype")
plot_genes_in_pseudotime(cds[keygenes,], color_by = "State")
plot_genes_in_pseudotime(cds[keygenes,], color_by = "time")
plot_genes_in_pseudotime(cds[keygenes,], color_by = "Pseudotime")

time_diff <- differentialGeneTest(cds[rownames(deg[1:2000,]),], fullModelFormulaStr = "~sm.ns(Pseudotime)")
time_diff = time_diff[time_diff$qval < 0.05,c(5,2,3,4,1,6,7)]

p = plot_pseudotime_heatmap(cds[time_diff$GeneName,], num_clusters = 5, cores = 1, 
                            show_rownames = F, return_heatmap = T)
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_heatmap.png', plot = p, height = 10, width = 4, bg = "transparent")
ggsave('../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_heatmap.svg', plot = p, height = 10, width = 4, bg = "transparent")

p$tree_row
clusters <- cutree(p$tree_row, k = 5)
clustering <- as.data.frame(clusters)
colnames(clustering) <- "Gene_Clusters"
clustering$gene = rownames(clustering)
table(clustering$Gene_Clusters)

# columns(org.Brapaoleracea.eg.db) 
# a = bitr(data1, 'GID', c("GENENAME", "GO", "GOALL", "ONTOLOGY", "ONTOLOGYALL"), org.Brapaoleracea.eg.db)

if (TRUE) {
  for (i in 1:5) {
    # i = 1
    data1 = rownames(clustering[clustering$Gene_Clusters == i,])
    data1[grepl('Bo', data1)] = paste0(data1[grepl('Bo', data1)], '.1')
    ego_ALL <- enrichGO(gene = data1, OrgDb = org.Brapaoleracea.eg.db, ont = "ALL", 
                        keyType = "GID", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
    write.table(as.data.frame(ego_ALL@result), paste("../revised_data/TT_8samples/monocle_res/TT_combined_Endo_5D7D8D9D11D_monocle_GO_c",i,".txt", sep = ''), quote = FALSE, row.names = TRUE, col.names = TRUE, sep="\t")
  }
  
  # 可视化
  if(TEUE){
    barplot(ego_ALL, x = "GeneRatio", color = "p.adjust", showCategory =10, split = "ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scales = "free", space = "free")
    
    dotplot(ego_ALL, x = "GeneRatio", color = "p.adjust", size = "Count", showCategory = 10, split="ONTOLOGY") + 
      facet_grid(ONTOLOGY~., scales = "free", space = "free") + NoLegend()
  }
}


