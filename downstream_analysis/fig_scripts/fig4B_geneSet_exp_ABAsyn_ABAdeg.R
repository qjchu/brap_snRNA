library(RColorBrewer)
library(reshape2)
library(Mfuzz)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggrepel)

sessionInfo()
setwd('scripts/')

rna_averages0 = read.table("../revised_data/CK3D/CK3D_RNA_harmony_clusters_res1_avg.txt")
colnames(rna_averages0) = paste('CK3D_',colnames(rna_averages0), sep = '')
rna_averages1 = read.table("../revised_data/CK5D_subcluster/CK5D_subcluster_RNA_harmony_clusters_res1_avg.txt")
colnames(rna_averages1) = paste('CK5D_',colnames(rna_averages1), sep = '')
rna_averages2 = read.table("../revised_data/CK7D/CK7D_RNA_harmony_clusters_res1_avg.txt")
colnames(rna_averages2) = paste('CK7D_',colnames(rna_averages2), sep = '')
rna_averages3 = read.table("../revised_data/CK9D/CK9D_RNA_harmony_clusters_res1_avg.txt")
colnames(rna_averages3) = paste('CK9D_',colnames(rna_averages3), sep = '')
rna_averages4 = read.table("../revised_data/CK11D_subcluster/CK11D_subcluster_RNA_harmony_clusters_res1_avg.txt")
colnames(rna_averages4) = paste('CK11D_',colnames(rna_averages4), sep = '')

rna_averages5 = read.table("../revised_data/TT5D/TT5D_combined_RNA_clusters_res1_Louvain_avg.txt")
colnames(rna_averages5) = paste('TT5D_',colnames(rna_averages5), sep = '')
rna_averages6 = read.table("../revised_data/TT7D/TT7D_combined_RNA_harmony_clusters_res1_5_avg.txt")
colnames(rna_averages6) = paste('TT7D_',colnames(rna_averages6), sep = '')
rna_averages7 = read.table("../revised_data/TT8D/TT8D_combined_RNA_harmony_clusters_res1_avg_rep1_rep3.txt")
colnames(rna_averages7) = paste('TT8D_',colnames(rna_averages7), sep = '')
rna_averages8 = read.table("../revised_data/TT9D/TT9D_combined_RNA_harmony_clusters_res1_Louvain_avg.txt")
colnames(rna_averages8) = paste('TT9D_',colnames(rna_averages8), sep = '')
rna_averages9 = read.table("../revised_data/TT11D/TT11D_combined_RNA_clusters_res1_avg.txt")
colnames(rna_averages9) = paste('TT11D_',colnames(rna_averages9), sep = '')

############################### function FindBchiHomo
FindBchiHomo = function(atha_gene, exp_data){
  # find homologous gene in bchi
  if(TRUE){
    transfer = read.table('bchi_to_atha_besthit.txt', sep = ' ')
    colnames(transfer) = c('bchi_id', 'atha_id')
    rownames(transfer) = transfer$bchi_id
    
    table(transfer$atha_id %in% atha_gene$atha_gene)
    bchi_gene = transfer[transfer$atha_id %in% atha_gene$atha_gene,]
    bchi_gene$symbol = atha_gene[bchi_gene$atha_id,2]
  }
  
  bchi_gene = bchi_gene[bchi_gene$bchi_id %in% rownames(exp_data),]
  data = exp_data[bchi_gene$bchi_id,]
  data = as.data.frame(data)
  data$gene = rownames(data)
  pdata <- reshape2::melt(data)
  colnames(pdata) = c('Gene', 'Time', 'Exp')
  pdata$sample = gsub('_.*','',pdata$Time)
  pdata$Label = c(rep(NA, length(pdata$Gene)-length(bchi_gene$symbol)), paste0(bchi_gene$bchi_id,'_',bchi_gene$symbol))
  
  return(pdata)
}

############################### ALL
if(TRUE){
  # CK ALL
  a = cbind(rna_averages0[rownames(rna_averages0),], 
            rna_averages1[rownames(rna_averages0),],
            rna_averages2[rownames(rna_averages0),],
            rna_averages3[rownames(rna_averages0),],
            rna_averages4[rownames(rna_averages0),])
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT ALL
  a = cbind(rna_averages5[rownames(rna_averages0),], 
            rna_averages6[rownames(rna_averages0),],
            rna_averages7[rownames(rna_averages0),],
            rna_averages8[rownames(rna_averages0),],
            rna_averages9[rownames(rna_averages0),])
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
}

############################### SC-SUS
if(TRUE){
  # CK SC-SUS
  a = cbind(rna_averages0[rownames(rna_averages0),16], 
            rna_averages1[rownames(rna_averages0),15],
            rna_averages2[rownames(rna_averages0),21])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT SC-SUS
  b = cbind(rna_averages5[rownames(rna_averages0),16], 
            rna_averages6[rownames(rna_averages0),24],
            rna_averages7[rownames(rna_averages0),23],
            rna_averages8[rownames(rna_averages0),26])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### PEN VS endosperm-pcd
if(TRUE){
  # CK PEN
  a = cbind(rna_averages0[rownames(rna_averages0),18], 
            rna_averages1[rownames(rna_averages0),3],
            rna_averages2[rownames(rna_averages0),3],
            rna_averages3[rownames(rna_averages0),9],
            rna_averages4[rownames(rna_averages0),6])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT PEN
  b = cbind(rna_averages5[rownames(rna_averages0),18], 
            rna_averages6[rownames(rna_averages0),11],
            rna_averages7[rownames(rna_averages0),15],
            rna_averages8[rownames(rna_averages0),20],
            rna_averages9[rownames(rna_averages0),17])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### PEN VS endosperm-active
if(TRUE){
  # CK PEN
  a = cbind(rna_averages0[rownames(rna_averages0),18], 
            rna_averages1[rownames(rna_averages0),3],
            rna_averages2[rownames(rna_averages0),3],
            rna_averages3[rownames(rna_averages0),9],
            rna_averages4[rownames(rna_averages0),6])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT endosperm-active
  b = cbind(rna_averages5[rownames(rna_averages0),24], 
            rna_averages6[rownames(rna_averages0),25],
            rna_averages7[rownames(rna_averages0),22],
            rna_averages8[rownames(rna_averages0),28])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### CZSC-phloem-xylem
if(TRUE){
  # CK CZSC-phloem-xylem
  a = cbind(rna_averages0[rownames(rna_averages0),19], 
            rna_averages1[rownames(rna_averages0),17],
            rna_averages2[rownames(rna_averages0),22],
            rna_averages3[rownames(rna_averages0),22],
            rna_averages4[rownames(rna_averages0),18])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT CZSC-phloem-xylem
  b = cbind(rna_averages5[rownames(rna_averages0),23], 
            rna_averages6[rownames(rna_averages0),23],
            rna_averages7[rownames(rna_averages0),21],
            rna_averages8[rownames(rna_averages0),22],
            rna_averages9[rownames(rna_averages0),13])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### SC-oi
if(TRUE){
  # CK SC-oi
  a = cbind(rna_averages0[rownames(rna_averages0),1], 
            rna_averages1[rownames(rna_averages0),1],
            rna_averages2[rownames(rna_averages0),1],
            rna_averages3[rownames(rna_averages0),12],
            rna_averages4[rownames(rna_averages0),14])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  # TT SC-oi
  b = cbind(rna_averages5[rownames(rna_averages0),3], 
            rna_averages6[rownames(rna_averages0),1],
            rna_averages7[rownames(rna_averages0),5],
            rna_averages8[rownames(rna_averages0),8],
            rna_averages9[rownames(rna_averages0),5])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### SC-ii
if(TRUE){
  a = cbind(rna_averages0[rownames(rna_averages0),7], 
            rna_averages1[rownames(rna_averages0),8],
            rna_averages2[rownames(rna_averages0),5],
            rna_averages3[rownames(rna_averages0),5],
            rna_averages4[rownames(rna_averages0),5])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  b = cbind(rna_averages5[rownames(rna_averages0),4],
            rna_averages6[rownames(rna_averages0),2],
            rna_averages7[rownames(rna_averages0),4],
            rna_averages8[rownames(rna_averages0),17],
            rna_averages9[rownames(rna_averages0),6])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### EP_Dividing_cell
if(TRUE){
  a = cbind(rna_averages0[rownames(rna_averages0),3], 
            rna_averages1[rownames(rna_averages0),11],
            rna_averages2[rownames(rna_averages0),11],
            rna_averages3[rownames(rna_averages0),10],
            rna_averages4[rownames(rna_averages0),13])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  b = cbind(rna_averages5[rownames(rna_averages0),15], 
            rna_averages6[rownames(rna_averages0),22],
            rna_averages7[rownames(rna_averages0),17],
            rna_averages8[rownames(rna_averages0),24])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

############################### CZSC
if(TRUE){
  a = cbind(rna_averages0[rownames(rna_averages0),5], 
            rna_averages1[rownames(rna_averages0),6],
            rna_averages2[rownames(rna_averages0),8],
            rna_averages3[rownames(rna_averages0),4],
            rna_averages4[rownames(rna_averages0),17])
  colnames(a) = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D')
  rownames(a) = rownames(rna_averages0)
  a = na.omit(a)
  
  b = cbind(rna_averages5[rownames(rna_averages0),8], 
            rna_averages6[rownames(rna_averages0),19],
            rna_averages7[rownames(rna_averages0),2],
            rna_averages8[rownames(rna_averages0),15],
            rna_averages8[rownames(rna_averages0),1])
  colnames(b) = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D')
  rownames(b) = rownames(rna_averages0)
  b = na.omit(b)
}

# ABA synthesis deepseeek
if(TRUE){
  atha_gene = c('AT5G67030','AT3G14440','AT1G30100','AT3G24220','AT1G78390','AT1G52340','AT1G16540','AT2G27150','AT1G60680','AT4G34000','AT3G19290')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('ABA1','NCED3','NCED5','NCED6','NCED9','ABA2','ABA3','AAO','LOS5/ABA4','ABF3','ABF4')
  rownames(atha_gene) = atha_gene$atha_gene
}

# ABA degradation deepseeek
if(TRUE){
  atha_gene = c('AT4G19230','AT2G29090','AT5G45340','AT3G19270','AT1G07260','AT1G78370','AT2G41040')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('CYP707A1','CYP707A2','CYP707A3','CYP707A4','UGT71B6','ABA-GT','DPA1')
  rownames(atha_gene) = atha_gene$atha_gene
}

# IAA synthesis deepseeek
if(TRUE){
  atha_gene = c('AT4G32540','AT4G13260','AT5G11320','AT1G70560','AT4G24670','AT4G39950','AT2G22330','AT3G44310','AT1G04610','AT5G43890','AT5G25620')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('YUC1','YUC2','YUC4','TAA1','TAR2','CYP79B2','CYP79B3','NIT1','YUC3','YUC5','YUC6')
  rownames(atha_gene) = atha_gene$atha_gene
}

# IAA degradation deepseeek
if(TRUE){
  atha_gene = c('AT2G14960','AT1G14130','AT1G14120','AT1G24100','AT1G24110')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('GH3.3','DAO1','DAO2','UGT74B1','UGT74D1')
  rownames(atha_gene) = atha_gene$atha_gene
}

# GA synthesis deepseeek
if(TRUE){
  atha_gene = c('AT4G02780','AT1G79460','AT5G25900','AT1G05160','AT5G51810','AT1G15550')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('GA1/CPS','GA2/KS','GA3/KO','GA4/KAO','GA20ox','GA3ox')
  rownames(atha_gene) = atha_gene$atha_gene
}

# GA degradation deepseeek
if(TRUE){
  atha_gene = c('AT1G78440','AT1G30040','AT1G02400','AT4G21200')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = c('GA2ox1','GA2ox2','GA2ox6','GA2ox8')
  rownames(atha_gene) = atha_gene$atha_gene
}

# XTH
if(TRUE){
  atha_gene = c('AT4G13080','AT4G13090','AT3G25050','AT2G06850','AT5G13870','AT5G65730','AT4G37800','AT1G11545','AT4G03210','AT2G14620','AT3G48580','AT5G57530','AT5G57540','AT4G25820','AT4G14130','AT3G23730','AT1G65310','AT4G30280','AT4G30290','AT5G48070','AT2G18800','AT5G57560','AT4G25810','AT4G30270','AT5G57550','AT4G28850','AT2G01850','AT1G14720','AT4G18990','AT1G32170','AT3G44990','AT2G36870','AT1G10550')
  atha_gene = as.data.frame(atha_gene)
  atha_gene$type = paste('XTH', seq(1,33), sep = '')
  rownames(atha_gene) = atha_gene$atha_gene
}

if(TRUE){
  pdata1 = FindBchiHomo(atha_gene = atha_gene, exp_data = a)
  pdata2 = FindBchiHomo(atha_gene = atha_gene, exp_data = b)
  wilcox.test(pdata1$Exp, pdata2$Exp, paired = FALSE)
  
  pdata = rbind(pdata1, pdata2)
  pdata$sample = substr(pdata$sample, 1, 2)
  ggplot(pdata, aes(x = Time, y = Exp, color = Time)) + 
    geom_boxplot(width = 0.5, cex = 1) + 
    geom_point(alpha = 0.6) +
    geom_line(aes(group = Gene), linetype=1, alpha = 0.6) +
    scale_color_manual(values = rep(c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22', '#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B'),50)) +
    facet_grid(~pdata$sample, space = 'free', scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill="white"),
          legend.position = "none") 
}

# ggsave CZSC
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSC_XTH.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave SC-SUS
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_GAsyn_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SUS_XTH.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave EP_Dividing_cell
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_EP_XTH.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave SCii
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_XTH.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_GAdegr_0.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCii_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  
}

# ggsave PEN VS endosperm-active
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_GAsyn_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoA_XTH.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave PEN VS endosperm-pcd
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_endoPCD_XTH_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave SCoi
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_XTH_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_IAAdegr_0.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_SCoi_ABAsyn_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}

# ggsave CZSC-phloem-xylem
if(TRUE){
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_XTH_00.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_GAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_GAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_IAAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_IAAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_ABAdegr.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
  ggsave('../revised_data/CK_TT_geneset/CK_TT_CZSCpx_ABAsyn.svg', plot = last_plot(), height = 3, width = 3.5, bg = "transparent")
}


# SC-oi endoA
if(TRUE){
  pdata1 = FindBchiHomo(atha_gene = atha_gene, exp_data = a)
  pdata2 = FindBchiHomo(atha_gene = atha_gene, exp_data = b)
  wilcox.test(pdata1$Exp, pdata2$Exp, paired = FALSE)
  
  pdata = rbind(pdata1, pdata2)
  pdata$sample = substr(pdata$sample, 1, 2)
  ggplot(pdata, aes(x = Time, y = Exp, color = Time)) + 
    geom_boxplot(width = 0.5, cex = 1) + 
    geom_point(alpha = 0.6) +
    geom_line(aes(group = Gene), linetype=1, alpha = 0.6) +
    scale_color_manual(values = rep(c( '#9ACD32', '#6B8E23','#556B2F','#32CD32','#228B22', '#BDB76B', '#FFD700','#FFA500','#DAA520','#B8860B'),50)) +
    facet_grid(~pdata$sample, space = 'free', scales = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill="white"),
          legend.position = "none") + 
    geom_text_repel(aes(label = Label), size = 3, direction = "y", hjust = 0, force = 10)
  
  # ggsave SCoi
  if(TRUE){
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_XTH_00.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_GAdegr.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_GAsyn.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_IAAdegr_0.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_IAAsyn.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_ABAdegr.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_SCoi_ABAsyn_00.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
  }
  
  # ggsave PEN VS endosperm-active
  if(TRUE){
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_ABAsyn.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_ABAdegr.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_IAAsyn.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_IAAdegr.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_GAsyn_00.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_GAdegr.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
    ggsave('../revised_data/CK_TT_geneset/labeled_CK_TT_endoA_XTH.svg', plot = last_plot(), height = 5, width = 15, bg = "transparent")
  }
  
}


