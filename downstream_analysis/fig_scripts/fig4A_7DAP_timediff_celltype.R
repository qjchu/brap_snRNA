library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(svglite)
library(ggpubr)
library(ggrepel)

setwd('D:/CHU THINKBOOK/A03 白菜修改/scripts/')

celltype_DEG = function(TT_data, CK_data, TT_celltype, CK_celltype, celltype){
  temp_data = cbind(TT_data[,TT_celltype], CK_data[rownames(TT_data),CK_celltype])
  colnames(temp_data) = c('TT', 'CK')
  rownames(temp_data) = rownames(TT_data)
  temp_data = as.data.frame(na.omit(temp_data))
  temp_data = temp_data[temp_data$TT + temp_data$CK > 0.5 & (temp_data$TT > 0) & (temp_data$CK > 0),]
  
  temp_data$logFC = log2(temp_data$TT/temp_data$CK)
  
  temp_data = temp_data[abs(temp_data$logFC)> 2,]
  temp_data$cluster = celltype
  temp_data$GeneID = rownames(temp_data)
  
  return(temp_data)
}

# TT7D VS CK7D 
if(TRUE){
  CK7D <- read.table('D:/CHU HP/A06 白菜单细胞/CK_batch3_10samples/CK7D/CK7D_RNA_harmony_clusters_res1_avg.txt', header = TRUE)
  colnames(CK7D) = c('CK7D SC-oi (c1)', 'CK7D SC (c2)', 'CK7D PEN (c3)', 'CK7D SC_Endosperm (c4)', 'CK7D SC-ii (c5)', 
                     'CK7D SC (c6)', 'CK7D SC_Dividing_Cell (c7)', 'CK7D CZSC (c8)', 'CK7D CZSC (c9)', 'CK7D SC (c10)', 
                     'CK7D EP_Dividing_Cell (c11)', 'CK7D SC (c12)','CK7D SC (c13)', 'CK7D SC-ii (c14)', 'CK7D SC_Endosperm (c15)', 
                     'CK7D CZE (c16)', 'CK7D EP_Dividing_Cell (c17)', 'CK7D SC-oi (c18)', 'CK7D MCE (c19)', 'CK7D SC (c20)', 
                     'CK7D SC_SUS (c21)', 'CK7D CZSC-phloem-xylem (c22)', 'CK7D SC (c23)', 'CK7D CZSC (c24)', 'CK7D EP (c25)')
  CK7D = CK7D[rowSums(CK7D) > 1,]
  
  TT7D <- read.table('../revised_data/TT7D/TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt', header = TRUE)
  colnames(TT7D) = c('TT7D SC-oi (c1)', 'TT7D SC-ii (c2)', 'TT7D SC (c3)', 'TT7D SC-endosperm (c4)', 'TT7D SC (c5)', 
                     'TT7D SC-ii (c6)', 'TT7D SC (c7)', 'TT7D SC (c8)', 'TT7D SC-dividing-cell (c9)', 'TT7D SC-oi (c10)', 
                     'TT7D Endosperm-PCD (c11)', 'TT7D SC (c12)', 'TT7D SC (c13)', 'TT7D SC-ii (c14)', 'TT7D SC (c15)', 
                     'TT7D CZSC (c16)', 'TT7D CZSC (c17)', 'TT7D SC (c18)', 'TT7D CZSC (c19)', 'TT7D SC-ii (c20)', 
                     'TT7D SC-oi (c21)', 'TT7D EP-dividing-cell (c22)', 'TT7D CZSC-phloem-xylem (c23)', 'TT7D SC-SUS (c24)', 'TT7D Endosperm-active (c25)')
  TT7D = TT7D[rowSums(TT7D) > 1,]
  
  data = cbind(TT7D, CK7D[rownames(TT7D),])
  dat = na.omit(data)
  
  # spearman corr
  a = sva::ComBat(scale(dat), batch = c(rep('TT7D',25), rep('CK7D',25)))
  M3 <- cor(a, method = 'spearman')
  pheatmap::pheatmap(M3[1:25,26:50])
  
  # 'TT7D Endosperm-active (c25)' 为正常功能的胚乳
  temp_data1 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D Endosperm-active (c25)', CK_celltype = 'CK7D PEN (c3)', celltype = 'Endosperm-active')
  temp_data2 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D CZSC (c16)', CK_celltype = 'CK7D CZSC (c9)', celltype = 'CZSC')
  temp_data3 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D CZSC-phloem-xylem (c23)', CK_celltype = 'CK7D CZSC-phloem-xylem (c22)', celltype = 'CZSC-phloem-xylem')
  temp_data4 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D EP-dividing-cell (c22)', CK_celltype = 'CK7D EP_Dividing_Cell (c17)', celltype = 'EP-dividing-cell')
  
  # TT7D Endosperm-PCD (c11) 为退化的胚乳，对应CK PEN进行比较
  temp_data5 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D Endosperm-PCD (c11)', CK_celltype = 'CK7D PEN (c3)', celltype = 'Endosperm-PCD')
  
  temp_data6 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC-oi (c1)', CK_celltype = 'CK7D SC-oi (c1)', celltype = 'SC-oi')
  temp_data7 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC-SUS (c24)', CK_celltype = 'CK7D SC_SUS (c21)', celltype = 'SC-SUS')
  temp_data8 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC-dividing-cell (c9)', CK_celltype = 'CK7D SC_Dividing_Cell (c7)', celltype = 'SC_Dividing_Cell')
  temp_data9 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC (c5)', CK_celltype = 'CK7D SC (c2)', celltype = 'SC1')
  temp_data10 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC (c7)', CK_celltype = 'CK7D SC (c10)', celltype = 'SC2')
  temp_data11 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC (c15)', CK_celltype = 'CK7D SC (c12)', celltype = 'SC3')
  temp_data12 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC-ii (c2)', CK_celltype = 'CK7D SC-ii (c5)', celltype = 'SC-ii')
  temp_data13 = celltype_DEG(TT_data = TT7D, CK_data = CK7D, TT_celltype = 'TT7D SC-endosperm (c4)', CK_celltype = 'CK7D SC_Endosperm (c15)', celltype = 'SC_Endosperm')
  
  # 所有比较结果
  if(TRUE){
    # 整合数据
    temp_data = rbind(temp_data1, temp_data2, temp_data3, temp_data4, temp_data5, temp_data6, temp_data7,temp_data8,temp_data9,temp_data10,temp_data11,temp_data12,temp_data13)
    
    temp_data$jittered_x <- jitter(as.numeric(factor(temp_data$cluster)), amount = 0.35)
    dfbar<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), 
                      y = c(max(temp_data1$logFC), max(temp_data2$logFC), max(temp_data3$logFC), max(temp_data4$logFC), max(temp_data5$logFC), max(temp_data6$logFC), max(temp_data7$logFC), max(temp_data8$logFC), max(temp_data9$logFC), max(temp_data10$logFC), max(temp_data11$logFC), max(temp_data12$logFC), max(temp_data13$logFC)))
    dfbar1<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), 
                       y = c(min(temp_data1$logFC), min(temp_data2$logFC), min(temp_data3$logFC), min(temp_data4$logFC), min(temp_data5$logFC), min(temp_data6$logFC), min(temp_data7$logFC), min(temp_data8$logFC), min(temp_data9$logFC), min(temp_data10$logFC), min(temp_data11$logFC), min(temp_data12$logFC), min(temp_data13$logFC)))
    dfcol<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), y = 0, 
                      label = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'))

    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC), size = 2) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = GeneID), force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x,y = y), height = 2, color = "white", fill = 'darkgreen', alpha = 0.6, show.legend = F) +
      geom_text(data = dfcol, aes(x = x, y=y, label=label), size = 5, color ="white") + 
      scale_color_manual(name=NULL, values = c("red","black"))+
      labs(x = "Comparison of celltype-specific DEGs between TT7D and CK7D", y = "log(FoldChange)") +
      theme_minimal() +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.line.y = element_line(color = "gray45", size = 1),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      )
    # ggsave('../CK_batch3_10samples/TT7D_VS_CK7D_celltype.svg', plot = last_plot(), height = 8, width = 24, bg = "transparent")
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC), size = 2,alpha = 0.6) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = GeneID), force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x,y = y), height = 3.4, color = "white", fill = 'darkgreen', alpha = 0.6, show.legend = F) +
      geom_text(data=dfcol, aes(x = x, y = y, label=label, angle = 90, vjust = 0.5, hjust = 0.5, size = 15), size = 5, color ="white") + 
      scale_color_manual(name=NULL, values = c("red","black"))+
      labs(# x = "Comparison of celltype-specific DEGs between TT7D and CK7D", 
        y = "log(FoldChange)") +
      theme_minimal() +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.line.x = element_line(color = "gray45", size = 1),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      ) + coord_flip()
    # ggsave('../CK_batch3_10samples/TT7D_VS_CK7D_celltype_v3.svg', plot = last_plot(), height = 24, width = 8, bg = "transparent")
    
    anno = read.table('../besthit_anno.txt', sep = '\t', fill = TRUE)
    colnames(anno) = c('bchi_gene', 'atha_gene', 'symbol', 'des1', 'des2')
    rownames(anno) = anno$bchi_gene
    anno$symbol[anno$symbol == 'None'] = anno$bchi_gene[anno$symbol == 'None']
    temp_data$symbol = anno[temp_data$GeneID,]$symbol
    temp_data$symbol[is.na(temp_data$symbol)] = temp_data$GeneID[is.na(temp_data$symbol)]
    temp_data$atha_gene = anno[temp_data$GeneID,]$atha_gene
    temp_data$des1 = anno[temp_data$GeneID,]$des1
    temp_data$des2 = anno[temp_data$GeneID,]$des2
    write.table(temp_data,file = '../CK_batch3_10samples/TT7D_VS_CK7D_celltype_DEG_v2.txt', sep = '\t', quote = FALSE)
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC), size = 2) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC,label = symbol), force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x,y = y), height = 2, color = "white", fill = mycol, alpha = 0.6, show.legend = F) +
      geom_text(data = dfcol, aes(x = x,y=y,label=label), size = 5, color ="white") + 
      scale_color_manual(name=NULL, values = c("red","black"))+
      labs(x = "Comparison of celltype-specific DEGs between TT7D and CK7D", y = "log(FoldChange)") +
      theme_minimal() +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.line.y = element_line(color = "gray45", size = 1),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      )
    # ggsave('../CK_batch3_10samples/TT7D_VS_CK7D_celltype_v2.svg', plot = last_plot(), height = 8, width = 24, bg = "transparent")
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC), size = 2,alpha = 0.6) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = symbol), force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x,y = y), height = 3.4, color = "white", fill = mycol, alpha = 0.6, show.legend = F) +
      geom_text(data=dfcol, aes(x = x, y = y, label=label, angle = 90, vjust = 0.5, hjust = 0.5, size = 15), size = 5, color ="white") + 
      scale_color_manual(name=NULL, values = c("red","black"))+
      labs(# x = "Comparison of celltype-specific DEGs between TT7D and CK7D", 
        y = "log(FoldChange)") +
      theme_minimal() +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.line.x = element_line(color = "gray45", size = 1),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      ) + coord_flip()
    # ggsave('../CK_batch3_10samples/TT7D_VS_CK7D_celltype_v4.svg', plot = last_plot(), height = 24, width = 8, bg = "transparent")
  }
  
  # fig4A 选取部分细胞类型进行展示
  if(TRUE){
    temp_data = rbind(temp_data1, temp_data2, temp_data3, temp_data4, temp_data5, temp_data6,temp_data7,temp_data12,temp_data13)
    temp_data$jittered_x <- jitter(as.numeric(factor(temp_data$cluster)), amount = 0.35)
    dfbar<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC-ii','SC_Endosperm'), 
                      y = c(max(temp_data1$logFC), max(temp_data2$logFC), max(temp_data3$logFC), max(temp_data4$logFC), max(temp_data5$logFC), max(temp_data6$logFC), max(temp_data7$logFC), max(temp_data12$logFC), max(temp_data13$logFC)))
    dfbar1<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC-ii','SC_Endosperm'), 
                       y = c(min(temp_data1$logFC), min(temp_data2$logFC), min(temp_data3$logFC), min(temp_data4$logFC), min(temp_data5$logFC), min(temp_data6$logFC), min(temp_data7$logFC), min(temp_data12$logFC), min(temp_data13$logFC)))
    dfcol<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC-ii','SC_Endosperm'), y = 0, 
                      label = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC-ii','SC_Endosperm'))
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC), size = 2,alpha = 0.6) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = GeneID), force = 1.2, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x,y = y), height = 3.4, color = "white", fill = 'darkgreen', alpha = 0.6, show.legend = F) +
      geom_text(data=dfcol, aes(x = x, y = y, label=label, angle = 90, vjust = 0.5, hjust = 0.5, size = 15), size = 5, color ="white") + 
      scale_color_manual(name=NULL, values = c("red","black"))+
      labs(# x = "Comparison of celltype-specific DEGs between TT7D and CK7D", 
        y = "Log2(FoldChange)") +
      theme_minimal() +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.line.x = element_line(color = "gray45", size = 1),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "top",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      ) + coord_flip()
    ggsave('../revised_data/TT7D_VS_CK7D_celltype_fig4A.svg', plot = last_plot(), height = 18, width = 8, bg = "transparent")
    
    anno = read.table('besthit_anno.txt', sep = '\t', fill = TRUE)
    colnames(anno) = c('bchi_gene', 'atha_gene', 'symbol', 'des1', 'des2')
    rownames(anno) = anno$bchi_gene
    anno$symbol[anno$symbol == 'None'] = anno$bchi_gene[anno$symbol == 'None']
    temp_data$symbol = anno[temp_data$GeneID,]$symbol
    temp_data$symbol[is.na(temp_data$symbol)] = temp_data$GeneID[is.na(temp_data$symbol)]
    temp_data$atha_gene = anno[temp_data$GeneID,]$atha_gene
    temp_data$des1 = anno[temp_data$GeneID,]$des1
    temp_data$des2 = anno[temp_data$GeneID,]$des2
    temp_data$highlight = NA
    temp_data$highlight[abs(temp_data$TT - temp_data$CK) > 1] = temp_data$logFC[abs(temp_data$TT - temp_data$CK) > 1]
    temp_data$symbol2 = NA
    # show = (temp_data$symbol %in% c('AGL62','ICE1','AGL91','UBP13','LEC1','TT2','CYP78A5','MYB65','MET1')) | (abs(temp_data$logFC) > 7) 
    show = temp_data$symbol %in% c('AGL62','ICE1','AGL91','UBP13','LEC1','TT2','CYP78A5','MYB65','MET1',
                                   'Bch06G010590','PAB8','XTH24','HB21','LCR84','DRM2','Bch04G013320','ENODL3','LCR30')
    temp_data$symbol2[show] = temp_data$symbol[show]
    
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC, color = highlight), size = 2, alpha = 0.96) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = symbol2), force = 3, box.padding = 0.8, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x, y = y), height = 3.6, color = "white", fill = 'darkgreen', alpha = 0.6, show.legend = F) +
      geom_text(data = dfcol, aes(x = x, y = y, label = label, angle = 90, vjust = 0.5, hjust = 0.5, size = 15), size = 5, color ="white") + 
      scale_color_gradient2(low = '#00e600', high = 'red')+
      labs(# x = "Comparison of celltype-specific DEGs between TT7D and CK7D", 
        y = "Log2(FoldChange)") +
      theme_minimal() + ylim(-14,18) +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.title.y = element_blank(),
            axis.line.x = element_line(color = "gray45", size = 1),
            axis.line.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      ) + coord_flip()
    ggsave('../revised_data/TT7D_VS_CK7D_celltype_fig4A_v2.svg', plot = last_plot(), height = 18, width = 8, bg = "transparent")
    ggsave('../revised_data/TT7D_VS_CK7D_celltype_fig4A_v2.png', plot = last_plot(), height = 18, width = 8, bg = "transparent")
    
    ## S8 展示所有细胞类型
    temp_data = rbind(temp_data1, temp_data2, temp_data3, temp_data4, temp_data5, temp_data6, temp_data7, temp_data8, temp_data9, temp_data10, temp_data11, temp_data12, temp_data13)
    
    temp_data$jittered_x <- jitter(as.numeric(factor(temp_data$cluster)), amount = 0.35)
    dfbar<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), 
                      y = c(max(temp_data1$logFC), max(temp_data2$logFC), max(temp_data3$logFC), max(temp_data4$logFC), max(temp_data5$logFC), max(temp_data6$logFC), max(temp_data7$logFC), max(temp_data8$logFC), max(temp_data9$logFC), max(temp_data10$logFC), max(temp_data11$logFC), max(temp_data12$logFC), max(temp_data13$logFC)))
    dfbar1<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), 
                       y = c(min(temp_data1$logFC), min(temp_data2$logFC), min(temp_data3$logFC), min(temp_data4$logFC), min(temp_data5$logFC), min(temp_data6$logFC), min(temp_data7$logFC), min(temp_data8$logFC), min(temp_data9$logFC), min(temp_data10$logFC), min(temp_data11$logFC), min(temp_data12$logFC), min(temp_data13$logFC)))
    dfcol<-data.frame(x = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'), y = 0, 
                      label = c('Endosperm-active','CZSC','CZSC-phloem-xylem','EP_Dividing_Cell','Endosperm-PCD','SC-oi','SC-SUS','SC_Dividing_Cell','SC1','SC2','SC3','SC-ii','SC_Endosperm'))
    
    anno = read.table('besthit_anno.txt', sep = '\t', fill = TRUE)
    colnames(anno) = c('bchi_gene', 'atha_gene', 'symbol', 'des1', 'des2')
    rownames(anno) = anno$bchi_gene
    anno$symbol[anno$symbol == 'None'] = anno$bchi_gene[anno$symbol == 'None']
    temp_data$symbol = anno[temp_data$GeneID,]$symbol
    temp_data$symbol[is.na(temp_data$symbol)] = temp_data$GeneID[is.na(temp_data$symbol)]
    temp_data$atha_gene = anno[temp_data$GeneID,]$atha_gene
    temp_data$des1 = anno[temp_data$GeneID,]$des1
    temp_data$des2 = anno[temp_data$GeneID,]$des2
    temp_data$highlight = NA
    temp_data$highlight[abs(temp_data$TT - temp_data$CK) > 1] = temp_data$logFC[abs(temp_data$TT - temp_data$CK) > 1]
    temp_data$symbol2 = NA
    show = (temp_data$symbol %in% c('AGL62','ICE1','AGL91','UBP13','LEC1','TT2','CYP78A5','MYB65','MET1')) | (abs(temp_data$logFC) > 5) 
    temp_data$symbol2[show] = temp_data$symbol[show]
    
    write.table(temp_data, file = '../revised_data/TT7D_VS_CK7D_celltype_tableS7.txt', sep = '\t', quote = FALSE, col.names = TRUE, row.names = TRUE)
    
    ggplot()+ 
      geom_col(data = dfbar, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_col(data = dfbar1, mapping = aes(x = x,y = y), fill = "#dcdcdc",alpha = 0.6) +
      geom_point(data = temp_data, aes(x = jittered_x, y = logFC, color = highlight), size = 2, alpha = 0.96) +
      geom_text_repel(data = temp_data, aes(x = jittered_x, y=logFC, label = symbol2), force = 3, box.padding = 0.3, arrow = arrow(length = unit(0.008, "npc"), type = "open", ends = "last")) +
      geom_tile(data = dfcol, aes(x = x, y = y), height = 3.6, color = "white", fill = 'darkgreen', alpha = 0.6, show.legend = F) +
      geom_text(data = dfcol, aes(x = x, y = y, label = label, angle = 0, vjust = 0.5, hjust = 0.5, size = 15), size = 5, color ="white") + 
      scale_color_gradient2(low = '#00e600', high = 'red')+
      labs(# x = "Comparison of celltype-specific DEGs between TT7D and CK7D", 
        y = "Log2(FoldChange)") +
      theme_minimal() + ylim(-14,18) +
      theme(axis.title = element_text(size = 13, color = "black"),
            axis.title.x = element_blank(),
            axis.line.y = element_line(color = "gray45", size = 1),
            axis.line.x = element_blank(),
            axis.text.x = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            legend.direction = "vertical",
            legend.justification = c(1,0),
            legend.text = element_text(size = 15)
      )
    ggsave('../revised_data/TT7D_VS_CK7D_celltype_fig4A_v3.svg', plot = last_plot(), height = 8, width = 24, bg = "transparent")
    ggsave('../revised_data/TT7D_VS_CK7D_celltype_fig4A_v3.png', plot = last_plot(), height = 8, width = 24, bg = "transparent")
    
  }
  
}
