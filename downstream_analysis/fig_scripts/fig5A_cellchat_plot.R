library(ggplot2)
library(ggsci)

getwd()
setwd('scripts/')

# CK3D
if(TRUE){
  data = read.table('../cell_to_cell_communicate/CK3D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  for (i in 1:20) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK3D'
  aa$cluster = paste('CK3D_', rownames(aa), sep = '')
  aa$celltype = c('SC-oi (c1)','SC-dividing-cell (c2)',"EP-dividing-cell (c3)","SC-dividing-cell (c4)","CZSC (c5)","SC (c6)","SC-ii (c7)","SC-dividing-cell (c8)","SC (c9)","SC-endosperm (c10)",
                  "SC-ii (c11)","SC-oi (c12)","SC (c13)", "SC-endosperm (c14)","SC (c15)","SC-SUS (c16)","SC (c17)","PEN/MCE (c18)","CZSC-phloem-xylem (c19)","CZE (c20)")
  aa$type = c('SC','SC','EP','SC','SC','SC','SC','SC','SC','SC','SC','SC','SC','SC','SC','SC','SC','Endosperm','SC','Endosperm')
}
plot_data = aa

# CK5D
if(TRUE){
  data = read.table('../cell_to_cell_communicate/CK5D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  for (i in 1:18) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK5D'
  aa$cluster = paste('CK5D_', rownames(aa), sep = '')
  aa$celltype = c("SC-oi (c1)", "SC-endosperm (c2)","PEN (c3)","SC-dividing-cell (c4)","SC-oi (c5)", "CZSC (c6)","SC (c7)", "SC-ii (c8)", "SC (c9)", 
                  "SC-ii (c10)", "EP-dividing-cell (c11)", "SC-endosperm (c12)", "SC (c13)", "CZE (c14)", "SC-SUS (c15)", "MCE (c16)", "CZSC-phloem-xylem (c17)", "SC (c18)")
  aa$type = c('SC','SC','Endosperm','SC','SC','SC','SC','SC','SC','SC','EP','SC','SC','Endosperm','SC','Endosperm','SC','SC')
}
plot_data = rbind(plot_data,aa)

# CK7D
if(TRUE){
  data = read.table('../cell_to_cell_communicate/CK7D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  for (i in 1:25) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK7D'
  aa$cluster = paste('CK7D_', rownames(aa), sep = '')
  aa$celltype = c("SC-oi (c1)", "SC (c2)", "PEN (c3)", "SC-endosperm (c4)", "SC-ii (c5)","SC-endosperm (c6)", "SC-dividing-cell (c7)", "CZSC (c8)", "CZSC (c9)", "SC (c10)", "EP-dividing-cell (c11)", "SC (c12)", "SC (c13)", "SC-ii (c14)",
                  "SC-endosperm (c15)", "CZE (c16)", "EP-dividing-cell (c17)", "SC-oi (c18)", "MCE (c19)", "SC (c20)", "SC-SUS (c21)","CZSC-phloem-xylem (c22)", "SC (c23)", "CZSC (c24)", "EP (c25)" )
  aa$type = c('SC','SC','Endosperm','SC','SC','SC','SC','SC','SC','SC','EP','SC','SC','SC','SC','Endosperm','EP','SC','Endosperm','SC','SC','SC','SC','SC','EP')
}
plot_data = rbind(plot_data,aa)
  
# CK9D
if(TRUE){
  data = read.table('../cell_to_cell_communicate/CK9D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  for (i in 1:23) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK9D'
  aa$cluster = paste('CK9D_', rownames(aa), sep = '')
  aa$celltype = c("SC-endosperm (c1)", "CZE (c2)", "SC (c3)", "CZSC (c4)", "SC-ii (c5)","SC-endosperm (c6)", "SC (c7)", "SC (c8)",
                  "PEN (c9)", "EP-dividing-cell (c10)", "PEN (c11)", "SC-oi (c12)", "SC-oi (c13)", "CZSC (c14)", "PEN (c15)", "PEN-around-EP (c16)",
                  "SC (c17)", "SC (c18)", "MCE-PEN (c19)", "PEN (c20)", "EP (c21)", "CZSC-phloem-xylem (c22)", "CZE (c23)")
  aa$type = c('SC','Endosperm','SC','SC','SC','SC','SC','SC','Endosperm','EP','Endosperm','SC','SC','SC','Endosperm','Endosperm','SC','SC','Endosperm','Endosperm','EP','SC','Endosperm')
}
plot_data = rbind(plot_data,aa)

# CK11D
if(TRUE){
  data = read.table('../cell_to_cell_communicate/CK11D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  for (i in 1:27) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK11D'
  aa$cluster = paste('CK11D_', rownames(aa), sep = '')
  aa$celltype = c("PEN (c1)", "PEN (c2)", "SC (c3)", "EP (c4)","SC-ii (c5)", "PEN (c6)", "EP (c7)","EP (c8)", "CZE (c9)", "SC (c10)",
                  "SC-PEN (c11)", "CZE (c12)", "EP-dividing-cell (c13)", "SC-oi (c14)", "EP (c15)","PEN (c16)", "CZSC (c17)", "CZSC-phloem-xylem (c18)", "CZE (c19)", "CZSC (c20)",
                  "EP (c21)", "PEN (c22)", "CZE (c23)","EP (c24)", "CZSC (c25)","SC-ii (c26)", "SC (c27)")
  aa$type = c('Endosperm','Endosperm','SC','EP','SC','Endosperm','EP','EP','Endosperm','SC','SC','Endosperm','EP','SC','EP','Endosperm','SC','SC','Endosperm','SC','EP','Endosperm','Endosperm','EP','SC','SC','SC')
}
plot_data = rbind(plot_data,aa)

plot_data_CK = plot_data
plot_data_CK$treat = 'CK'

plot_data$time = factor(plot_data$time, levels = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D'))
ggplot(plot_data, aes(x = aa_source, y = aa_target, size = aa_self, color = type)) +
  geom_point() +
  geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 1.8) +
  scale_color_npg("nrc") + 
  scale_size_continuous(range = c(2,15)) +
  xlim(0,350) + ylim(0,350) +
  facet_grid(~plot_data$time, scales = 'free', space = 'free') +
  theme(text = element_text(size = 15),
        axis.line = element_line(color = 'darkgray'),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill="white"),
        legend.position = "right")
ggsave('../CK_batch3_10samples/CK_cellchat.svg', plot = last_plot(), height = 4.3, width = 20)
ggsave('../CK_batch3_10samples/CK_cellchat.png', plot = last_plot(), height = 4.3, width = 20, bg = "transparent")

# ggplot(plot_data, aes(x = aa_source, y = aa_target, size = aa_self, color = time)) +
#   geom_point() +
#   geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 3) +
#   # scale_color_igv("default") + 
#   scale_color_manual(values = c('#9ACD32','#6B8E23','#556B2F','#32CD32','#228B22')) +
#   scale_size_continuous(range = c(1,10)) +
#   xlim(0,350) + ylim(0,350) +
#   theme(text = element_text(size = 15),
#         axis.line = element_line(color = 'darkgray'),
#         panel.background = element_rect(fill = 'white'),
#         plot.background = element_rect(fill="white"),
#         legend.position = "right")


################ TT
# TT5D
if(TRUE){
  data = read.table('../revised_data/TT5D/TT5D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  
  table(data$source)
  for (i in 1:24) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK3D'
  aa$cluster = paste('TT5D_', rownames(aa), sep = '')
  aa$celltype = c('SC (c1)', 'SC-endosperm (c2)', 'SC-oi (c3)', 'SC-ii (c4)', 'SC-dividing-cell (c5)', 'SC-endosperm (c6)', 'SC (c7)', 'CZSC (c8)', 'SC (c9)', 'CZSC (c10)', 'SC-dividing-cell (c11)', 'SC (c12)', 'SC-ii (c13)', 'CZSC (c14)', 'EP-dividing-cell (c15)', 'SC-SUS (c16)', 'CZSC (c17)', 'Endosperm-PCD (c18)', 'SC-dividing-cell (c19)', 'SC-dividing-cell (c20)', 'SC (c21)', 'SC (c22)', 'CZSC-phloem-xylem (c23)', 'Endosperm-active (c24)')
  aa$type = c(rep('SC',14),'EP','SC','SC','Endosperm',rep('SC',5),'Endosperm')
}
plot_data = aa

# TT7D
if(TRUE){
  data = read.table('../revised_data/TT7D/TT7D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  
  table(data$source)
  for (i in 1:25) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK5D'
  aa$cluster = paste('TT7D_', rownames(aa), sep = '')
  aa$celltype = c('SC-oi (c1)', 'SC-ii (c2)', 'SC (c3)', 'SC-endosperm (c4)', 'SC (c5)', 'SC-ii (c6)', 'SC (c7)', 'SC (c8)', 'SC-dividing-cell (c9)', 'SC-oi (c10)', 'Endosperm-PCD (c11)', 'SC (c12)', 'SC (c13)', 'SC-ii (c14)', 'SC (c15)', 'CZSC (c16)', 'CZSC (c17)', 'SC (c18)', 'CZSC (c19)', 'SC-ii (c20)', 'SC-oi (c21)', 'EP-dividing-cell (c22)', 'CZSC-phloem-xylem (c23)', 'SC-SUS (c24)', 'Endosperm-active (c25)')
  aa$type = c(rep('SC',10),'Endosperm',rep('SC',10),'EP','SC','SC','Endosperm')
}
plot_data = rbind(plot_data,aa)

# TT8D
if(TRUE){
  data = read.table('../revised_data/TT8D/TT8D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  
  table(data$source)
  for (i in 1:24) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK7D'
  aa$cluster = paste('TT8D_', rownames(aa), sep = '')
  aa$celltype = c('SC (c1)', 'CZSC (c2)', 'SC (c3)', 'SC-ii (c4)', 'SC-oi (c5)', 'SC-endosperm (c6)', 'SC-dividing-cell (c7)', 'SC (c8)', 'CZSC (c9)', 'SC-oi (c10)', 'SC (c11)', 'CZSC (c12)', 'SC-ii (c13)', 'SC (c14)', 'Endosperm-PCD (c15)', 'SC (c16)', 'EP-dividing-cell (c17)', 'SC (c18)', 'Endosperm-PCD (c19)', 'SC (c20)', 'CZSC-phloem-xylem (c21)', 'Endosperm-active (c22)', 'SC-SUS (c23)', 'Endosperm-active (c24)')
  aa$type = c(rep('SC',14),'Endosperm','SC','EP','SC','Endosperm','SC','SC','Endosperm','SC','Endosperm')
}
plot_data = rbind(plot_data,aa)

# TT9D
if(TRUE){
  data = read.table('../revised_data/TT9D/TT9D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  
  table(data$source)
  for (i in 1:28) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK9D'
  aa$cluster = paste('TT9D_', rownames(aa), sep = '')
  aa$celltype = c('SC (c1)', 'SC (c2)', 'SC (c3)', 'SC-endosperm (c4)', 'CZSC (c5)', 'SC (c6)', 'SC-ii (c7)', 'SC-oi (c8)', 'SC-oi (c9)', 'SC-ii (c10)', 'SC (c11)', 'CZSC (c12)', 'CZSC (c13)', 'SC (c14)', 'CZSC (c15)', 'SC (c16)', 'SC (c17)', 'SC-ii (c18)', 'SC-dividing-cell (c19)', 'Endosperm-PCD (c20)', 'CZSC (c21)', 'CZSC-phloem-xylem (c22)', 'SC (c23)', 'EP-dividing-cell (c24)', 'SC (c25)', 'SC-SUS (c26)', 'SC (c27)', 'Endosperm-active (c28)')
  aa$type = c(rep('SC',19),'Endosperm','SC','SC','SC','EP','SC','SC','SC','Endosperm')
}
plot_data = rbind(plot_data,aa)

# TT11D
if(TRUE){
  data = read.table('../revised_data/TT11D/TT11D_net_lr.txt', header = TRUE, sep = ',')
  aa_source = c()
  aa_target = c()
  aa_self = c()
  
  table(data$source)
  for (i in 1:17) {
    aa_self = c(aa_self, length(data[data$target == i & data$source == i,1]))
    aa_source = c(aa_source, length(data[data$source == i,1])-length(data[data$target == i & data$source == i,1]))
    aa_target = c(aa_target, length(data[data$target == i,1])-length(data[data$target == i & data$source == i,1]))
  }
  
  aa = cbind(aa_self, aa_source, aa_target)
  aa = as.data.frame(aa)
  aa$time = 'CK11D'
  aa$cluster = paste('TT11D_', rownames(aa), sep = '')
  aa$celltype = c('CZSC (c1)', 'SC (c2)', 'SC (c3)', 'SC (c4)', 'SC-oi (c5)', 'SC (c6)', 'SC (c7)', 'SC (c8)', 'CZSC (c9)', 'SC (c10)', 'CZSC (c11)', 'SC (c12)', 'CZSC-phloem-xylem (c13)', 'SC (c14)', 'SC (c15)', 'SC (c16)', 'Endosperm (c17)')
  aa$type = c(rep('SC',16), 'Endosperm')
}
plot_data = rbind(plot_data,aa)

plot_data_TT = plot_data
plot_data_TT$treat = 'TT'

plot_data$time = factor(plot_data$time, levels = c('TT5D', 'TT7D', 'TT8D', 'TT9D', 'TT11D'))
ggplot(plot_data, aes(x = aa_source, y = aa_target, size = aa_self, color = type)) +
  geom_point() +
  geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 1.8) +
  scale_color_npg("nrc") + 
  # scale_color_manual(values = c('#9ACD32','#6B8E23','#556B2F','#32CD32','#228B22')) +
  scale_size_continuous(range = c(2,12)) +
  xlim(0,350) + ylim(0,350) +
  facet_grid(~plot_data$time, scales = 'free', space = 'free') +
  theme(text = element_text(size = 15),
        axis.line = element_line(color = 'darkgray'),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill="white"),
        legend.position = "right")
ggsave('../TT_batch3_8samples/TT_cellchat.svg', plot = last_plot(), height = 4.3, width = 20)
ggsave('../TT_batch3_8samples/TT_cellchat.png', plot = last_plot(), height = 4.3, width = 20, bg = "transparent")

## CK_TT
plot_data_all = rbind(plot_data_CK, plot_data_TT)
plot_data_all$time = factor(plot_data_all$time, levels = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D'))
ggplot(plot_data_all, aes(x = aa_source, y = aa_target, size = aa_self, color = type)) +
  geom_point() +
  geom_text(aes(label = celltype), vjust = 0.5, hjust = 0.5, color = "gray20", size = 1.8) +
  scale_color_npg("nrc") + 
  scale_size_continuous(range = c(2,15)) +
  xlim(0,350) + ylim(0,350) +
  facet_grid(plot_data_all$treat~plot_data_all$time, scales = 'free', space = 'free') +
  theme(text = element_text(size = 15),
        axis.line = element_line(color = 'darkgray'),
        panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill="white"),
        legend.position = "right")
ggsave('../revised_data/CK_TT_cellchat.svg', plot = last_plot(), height = 8.1, width = 20)
ggsave('../revised_data/CK_TT_cellchat.png', plot = last_plot(), height = 8.1, width = 20, bg = "transparent")


