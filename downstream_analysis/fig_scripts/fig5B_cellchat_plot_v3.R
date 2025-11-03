library(CellChat)
library(Seurat)
library(tidyverse)
library(NMF)
library(ggalluvial)
library(patchwork)
library(svglite)
library(ggrepel)
library(patchwork)

library(ggplot2)
library(ggsci)
library(dplyr)

getwd()
setwd('D:/CHU THINKBOOK/A03 白菜修改/scripts/')

if(TRUE){
  # CK3D
  if(TRUE){
    data = read.table('../revised_data/CK3D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                  "1" = 'SC-oi (c1)',
                  "2" = 'SC-dividing-cell (c2)',
                  "3" = "EP-dividing-cell (c3)",
                  "4" = "SC-dividing-cell (c4)",
                  "5" = "CZSC (c5)",
                  "6" = "SC (c6)",
                  "7" = "SC-ii (c7)",
                  "8" = "SC-dividing-cell (c8)",
                  "9" = "SC (c9)",
                  "10" = "SC-endosperm (c10)",
                  "11" = "SC-ii (c11)",
                  "12" = "SC-oi (c12)",
                  "13" = "SC (c13)",
                  "14" = "SC-endosperm (c14)",
                  "15" = "SC (c15)",
                  "16" = "SC-SUS (c16)",
                  "17" = "SC (c17)",
                  "18" = "PEN/MCE (c18)",
                  "19" = "CZSC-phloem-xylem (c19)",
                  "20" = "CZE (c20)")
    data$target2 <- recode(data$target, 
                           "1" = 'SC-oi (c1)',
                           "2" = 'SC-dividing-cell (c2)',
                           "3" = "EP-dividing-cell (c3)",
                           "4" = "SC-dividing-cell (c4)",
                           "5" = "CZSC (c5)",
                           "6" = "SC (c6)",
                           "7" = "SC-ii (c7)",
                           "8" = "SC-dividing-cell (c8)",
                           "9" = "SC (c9)",
                           "10" = "SC-endosperm (c10)",
                           "11" = "SC-ii (c11)",
                           "12" = "SC-oi (c12)",
                           "13" = "SC (c13)",
                           "14" = "SC-endosperm (c14)",
                           "15" = "SC (c15)",
                           "16" = "SC-SUS (c16)",
                           "17" = "SC (c17)",
                           "18" = "PEN/MCE (c18)",
                           "19" = "CZSC-phloem-xylem (c19)",
                           "20" = "CZE (c20)")
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = gsub(".*PEN.*", "Endosperm", data$source3)
    data$target4 = gsub(".*PEN.*", "Endosperm", data$target3)
    data$source4 = gsub(".*CZE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*CZE.*", "Endosperm", data$target4)
    data$source4 = gsub(".*MCE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*MCE.*", "Endosperm", data$target4)
    # plot_data1 = data[data$interaction_name %in% c('AT1G62340->AT5G44700','AT4G15800->AT2G46330','AT4G15800->AT3G61640','AT2G38540->AT5G37780','CYSB->MSS1'),]
    plot_data1 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data1$time = 'CK3D'
  }
  
  # CK5D
  if(TRUE){
    data = read.table('../revised_data/CK5D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = "SC-oi (c1)",
                           "2" = "SC-endosperm (c2)",
                           "3" = "PEN (c3)",
                           "4" = "SC-dividing-cell (c4)",
                           "5" = "SC-oi (c5)",
                           "6" = "CZSC (c6)",
                           "7" = "SC (c7)",
                           "8" = "SC-ii (c8)",
                           "9" = "SC (c9)",
                           "10" = "SC-ii (c10)",
                           "11" = "EP-dividing-cell (c11)",
                           "12" = "SC-endosperm (c12)",
                           "13" = "SC (c13)",
                           "14" = "CZE (c14)",
                           "15" = "SC-SUS (c15)",
                           "16" = "MCE (c16)",
                           "17" = "CZSC-phloem-xylem (c17)",
                           "18" = "SC (c18)")
    data$target2 <- recode(data$target, 
                           "1" = "SC-oi (c1)",
                           "2" = "SC-endosperm (c2)",
                           "3" = "PEN (c3)",
                           "4" = "SC-dividing-cell (c4)",
                           "5" = "SC-oi (c5)",
                           "6" = "CZSC (c6)",
                           "7" = "SC (c7)",
                           "8" = "SC-ii (c8)",
                           "9" = "SC (c9)",
                           "10" = "SC-ii (c10)",
                           "11" = "EP-dividing-cell (c11)",
                           "12" = "SC-endosperm (c12)",
                           "13" = "SC (c13)",
                           "14" = "CZE (c14)",
                           "15" = "SC-SUS (c15)",
                           "16" = "MCE (c16)",
                           "17" = "CZSC-phloem-xylem (c17)",
                           "18" = "SC (c18)")
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = gsub(".*PEN.*", "Endosperm", data$source3)
    data$target4 = gsub(".*PEN.*", "Endosperm", data$target3)
    data$source4 = gsub(".*CZE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*CZE.*", "Endosperm", data$target4)
    data$source4 = gsub(".*MCE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*MCE.*", "Endosperm", data$target4)
    plot_data2 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data2$time = 'CK5D'
  }
  
  # CK7D
  if(TRUE){
    data = read.table('../revised_data/CK7D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = "SC-oi (c1)",
                           "2" = "SC (c2)",
                           "3" = "PEN (c3)",
                           "4" = "SC-endosperm (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "SC-endosperm (c6)",
                           "7" = "SC-dividing-cell (c7)",
                           "8" = "CZSC (c8)",
                           "9" = "CZSC (c9)",
                           "10" = "SC (c10)",
                           "11" = "EP-dividing-cell (c11)",
                           "12" = "SC (c12)",
                           "13" = "SC (c13)",
                           "14" = "SC-ii (c14)",
                           "15" = "SC-endosperm (c15)",
                           "16" = "CZE (c16)",
                           "17" = "EP-dividing-cell (c17)",
                           "18" = "SC-oi (c18)",
                           "19" = "MCE (c19)",
                           "20" = "SC (c20)",
                           "21" = "SC-SUS (c21)",
                           "22" = "CZSC-phloem-xylem (c22)",
                           "23" = "SC (c23)",
                           "24" = "CZSC (c24)",
                           "25" = "EP (c25)")
    data$target2 <- recode(data$target, 
                           "1" = "SC-oi (c1)",
                           "2" = "SC (c2)",
                           "3" = "PEN (c3)",
                           "4" = "SC-endosperm (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "SC-endosperm (c6)",
                           "7" = "SC-dividing-cell (c7)",
                           "8" = "CZSC (c8)",
                           "9" = "CZSC (c9)",
                           "10" = "SC (c10)",
                           "11" = "EP-dividing-cell (c11)",
                           "12" = "SC (c12)",
                           "13" = "SC (c13)",
                           "14" = "SC-ii (c14)",
                           "15" = "SC-endosperm (c15)",
                           "16" = "CZE (c16)",
                           "17" = "EP-dividing-cell (c17)",
                           "18" = "SC-oi (c18)",
                           "19" = "MCE (c19)",
                           "20" = "SC (c20)",
                           "21" = "SC-SUS (c21)",
                           "22" = "CZSC-phloem-xylem (c22)",
                           "23" = "SC (c23)",
                           "24" = "CZSC (c24)",
                           "25" = "EP (c25)")
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = gsub(".*PEN.*", "Endosperm", data$source3)
    data$target4 = gsub(".*PEN.*", "Endosperm", data$target3)
    data$source4 = gsub(".*CZE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*CZE.*", "Endosperm", data$target4)
    data$source4 = gsub(".*MCE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*MCE.*", "Endosperm", data$target4)
    plot_data3 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data3$time = 'CK7D'
  }
  
  # CK9D
  if(TRUE){
    data = read.table('../revised_data/CK9D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = "SC-endosperm (c1)",
                           "2" = "CZE (c2)",
                           "3" = "SC (c3)",
                           "4" = "CZSC (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "SC-endosperm (c6)",
                           "7" = "SC (c7)",
                           "8" = "SC (c8)",
                           "9" = "PEN (c9)",
                           "10" = "EP-dividing-cell (c10)",
                           "11" = "PEN (c11)",
                           "12" = "SC-oi (c12)",
                           "13" = "SC-oi (c13)",
                           "14" = "CZSC (c14)",
                           "15" = "PEN (c15)",
                           "16" = "PEN-around-EP (c16)",
                           "17" = "SC (c17)",
                           "18" = "SC (c18)",
                           "19" = "MCE-PEN (c19)",
                           "20" = "PEN (c20)",
                           "21" = "EP (c21)",
                           "22" = "CZSC-phloem-xylem (c22)",
                           "23" = "CZE (c23)")
    data$target2 <- recode(data$target, 
                           "1" = "SC-endosperm (c1)",
                           "2" = "CZE (c2)",
                           "3" = "SC (c3)",
                           "4" = "CZSC (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "SC-endosperm (c6)",
                           "7" = "SC (c7)",
                           "8" = "SC (c8)",
                           "9" = "PEN (c9)",
                           "10" = "EP-dividing-cell (c10)",
                           "11" = "PEN (c11)",
                           "12" = "SC-oi (c12)",
                           "13" = "SC-oi (c13)",
                           "14" = "CZSC (c14)",
                           "15" = "PEN (c15)",
                           "16" = "PEN-around-EP (c16)",
                           "17" = "SC (c17)",
                           "18" = "SC (c18)",
                           "19" = "MCE-PEN (c19)",
                           "20" = "PEN (c20)",
                           "21" = "EP (c21)",
                           "22" = "CZSC-phloem-xylem (c22)",
                           "23" = "CZE (c23)")
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = gsub(".*PEN.*", "Endosperm", data$source3)
    data$target4 = gsub(".*PEN.*", "Endosperm", data$target3)
    data$source4 = gsub(".*CZE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*CZE.*", "Endosperm", data$target4)
    data$source4 = gsub(".*MCE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*MCE.*", "Endosperm", data$target4)
    plot_data4 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data4$time = 'CK9D'
  }
  
  # CK11D
  if(TRUE){
    data = read.table('../revised_data/CK11D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = "PEN (c1)",
                           "2" = "PEN (c2)",
                           "3" = "SC (c3)",
                           "4" = "EP (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "PEN (c6)",
                           "7" = "EP (c7)",
                           "8" = "EP (c8)",
                           "9" = "CZE (c9)",
                           "10" = "SC (c10)",
                           "11" = "SC-PEN (c11)",
                           "12" = "CZE (c12)",
                           "13" = "EP-dividing-cell (c13)",
                           "14" = "SC-oi (c14)",
                           "15" = "EP (c15)",
                           "16" = "PEN (c16)",
                           "17" = "CZSC (c17)",
                           "18" = "CZSC-phloem-xylem (c18)",
                           "19" = "CZE (c19)",
                           "20" = "CZSC (c20)",
                           "21" = "EP (c21)",
                           "22" = "PEN (c22)",
                           "23" = "CZE (c23)",
                           "24" = "EP (c24)",
                           "25" = "CZSC (c25)",
                           "26" = "SC-ii (c26)",
                           "27" = "SC (c27)")
    data$target2 <- recode(data$target, 
                           "1" = "PEN (c1)",
                           "2" = "PEN (c2)",
                           "3" = "SC (c3)",
                           "4" = "EP (c4)",
                           "5" = "SC-ii (c5)",
                           "6" = "PEN (c6)",
                           "7" = "EP (c7)",
                           "8" = "EP (c8)",
                           "9" = "CZE (c9)",
                           "10" = "SC (c10)",
                           "11" = "SC-PEN (c11)",
                           "12" = "CZE (c12)",
                           "13" = "EP-dividing-cell (c13)",
                           "14" = "SC-oi (c14)",
                           "15" = "EP (c15)",
                           "16" = "PEN (c16)",
                           "17" = "CZSC (c17)",
                           "18" = "CZSC-phloem-xylem (c18)",
                           "19" = "CZE (c19)",
                           "20" = "CZSC (c20)",
                           "21" = "EP (c21)",
                           "22" = "PEN (c22)",
                           "23" = "CZE (c23)",
                           "24" = "EP (c24)",
                           "25" = "CZSC (c25)",
                           "26" = "SC-ii (c26)",
                           "27" = "SC (c27)")
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = gsub(".*PEN.*", "Endosperm", data$source3)
    data$target4 = gsub(".*PEN.*", "Endosperm", data$target3)
    data$source4 = gsub(".*CZE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*CZE.*", "Endosperm", data$target4)
    data$source4 = gsub(".*MCE.*", "Endosperm", data$source4)
    data$target4 = gsub(".*MCE.*", "Endosperm", data$target4)
    plot_data5 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data5$time = 'CK11D'
  }
  
  # TT5D
  if(TRUE){
    data = read.table('../revised_data/TT5D/TT5D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = 'SC (c1)',
                           "2" = 'SC-endosperm (c2)',
                           "3" = 'SC-oi (c3)',
                           "4" = 'SC-ii (c4)',
                           "5" = 'SC-dividing-cell (c5)',
                           "6" = 'SC-endosperm (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'CZSC (c8)',
                           "9" = 'SC (c9)',
                           "10" = 'CZSC (c10)',
                           "11" = 'SC-dividing-cell (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'SC-ii (c13)',
                           "14" = 'CZSC (c14)',
                           "15" = 'EP-dividing-cell (c15)',
                           "16" = 'SC-SUS (c16)',
                           "17" = 'CZSC (c17)',
                           "18" = 'Endosperm-PCD (c18)',
                           "19" = 'SC-dividing-cell (c19)',
                           "20" = 'SC-dividing-cell (c20)',
                           "21" = 'SC (c21)',
                           "22" = 'SC (c22)',
                           "23" = 'CZSC-phloem-xylem (c23)',
                           "24" = 'Endosperm-active (c24)')
    data$target2 <- recode(data$target, 
                           "1" = 'SC (c1)',
                           "2" = 'SC-endosperm (c2)',
                           "3" = 'SC-oi (c3)',
                           "4" = 'SC-ii (c4)',
                           "5" = 'SC-dividing-cell (c5)',
                           "6" = 'SC-endosperm (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'CZSC (c8)',
                           "9" = 'SC (c9)',
                           "10" = 'CZSC (c10)',
                           "11" = 'SC-dividing-cell (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'SC-ii (c13)',
                           "14" = 'CZSC (c14)',
                           "15" = 'EP-dividing-cell (c15)',
                           "16" = 'SC-SUS (c16)',
                           "17" = 'CZSC (c17)',
                           "18" = 'Endosperm-PCD (c18)',
                           "19" = 'SC-dividing-cell (c19)',
                           "20" = 'SC-dividing-cell (c20)',
                           "21" = 'SC (c21)',
                           "22" = 'SC (c22)',
                           "23" = 'CZSC-phloem-xylem (c23)',
                           "24" = 'Endosperm-active (c24)')
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = sub("-PCD", "", data$source3)
    data$target4 = sub("-PCD", "", data$target3)
    data$source4 = sub("-active", "", data$source4)
    data$target4 = sub("-active", "", data$target4)
    plot_data6 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data6$time = 'TT5D'
  }
  
  # TT7D
  if(TRUE){
    data = read.table('../revised_data/TT7D/TT7D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = 'SC-oi (c1)',
                           "2" = 'SC-ii (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-endosperm (c4)',
                           "5" = 'SC (c5)',
                           "6" = 'SC-ii (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'SC-dividing-cell (c9)',
                           "10" = 'SC-oi (c10)',
                           "11" = 'Endosperm-PCD (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'SC (c13)',
                           "14" = 'SC-ii (c14)',
                           "15" = 'SC (c15)',
                           "16" = 'CZSC (c16)',
                           "17" = 'CZSC (c17)',
                           "18" = 'SC (c18)',
                           "19" = 'CZSC (c19)',
                           "20" = 'SC-ii (c20)',
                           "21" = 'SC-oi (c21)',
                           "22" = 'EP-dividing-cell (c22)',
                           "23" = 'CZSC-phloem-xylem (c23)',
                           "24" = 'SC-SUS (c24)',
                           "25" = 'Endosperm-active (c25)')
    data$target2 <- recode(data$target, 
                           "1" = 'SC-oi (c1)',
                           "2" = 'SC-ii (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-endosperm (c4)',
                           "5" = 'SC (c5)',
                           "6" = 'SC-ii (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'SC-dividing-cell (c9)',
                           "10" = 'SC-oi (c10)',
                           "11" = 'Endosperm-PCD (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'SC (c13)',
                           "14" = 'SC-ii (c14)',
                           "15" = 'SC (c15)',
                           "16" = 'CZSC (c16)',
                           "17" = 'CZSC (c17)',
                           "18" = 'SC (c18)',
                           "19" = 'CZSC (c19)',
                           "20" = 'SC-ii (c20)',
                           "21" = 'SC-oi (c21)',
                           "22" = 'EP-dividing-cell (c22)',
                           "23" = 'CZSC-phloem-xylem (c23)',
                           "24" = 'SC-SUS (c24)',
                           "25" = 'Endosperm-active (c25)')
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = sub("-PCD", "", data$source3)
    data$target4 = sub("-PCD", "", data$target3)
    data$source4 = sub("-active", "", data$source4)
    data$target4 = sub("-active", "", data$target4)
    plot_data7 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data7$time = 'TT7D'
  }
  
  # TT8D
  if(TRUE){
    data = read.table('../revised_data/TT8D/TT8D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = 'SC (c1)',
                           "2" = 'CZSC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-ii (c4)',
                           "5" = 'SC-oi (c5)',
                           "6" = 'SC-endosperm (c6)',
                           "7" = 'SC-dividing-cell (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'CZSC (c9)',
                           "10" = 'SC-oi (c10)',
                           "11" = 'SC (c11)',
                           "12" = 'CZSC (c12)',
                           "13" = 'SC-ii (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'Endosperm-PCD (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'EP-dividing-cell (c17)',
                           "18" = 'SC (c18)',
                           "19" = 'Endosperm-PCD (c19)',
                           "20" = 'SC (c20)',
                           "21" = 'CZSC-phloem-xylem (c21)',
                           "22" = 'Endosperm-active (c22)',
                           "23" = 'SC-SUS (c23)',
                           "24" = 'Endosperm-active (c24)')
    data$target2 <- recode(data$target, 
                           "1" = 'SC (c1)',
                           "2" = 'CZSC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-ii (c4)',
                           "5" = 'SC-oi (c5)',
                           "6" = 'SC-endosperm (c6)',
                           "7" = 'SC-dividing-cell (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'CZSC (c9)',
                           "10" = 'SC-oi (c10)',
                           "11" = 'SC (c11)',
                           "12" = 'CZSC (c12)',
                           "13" = 'SC-ii (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'Endosperm-PCD (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'EP-dividing-cell (c17)',
                           "18" = 'SC (c18)',
                           "19" = 'Endosperm-PCD (c19)',
                           "20" = 'SC (c20)',
                           "21" = 'CZSC-phloem-xylem (c21)',
                           "22" = 'Endosperm-active (c22)',
                           "23" = 'SC-SUS (c23)',
                           "24" = 'Endosperm-active (c24)')
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = sub("-PCD", "", data$source3)
    data$target4 = sub("-PCD", "", data$target3)
    data$source4 = sub("-active", "", data$source4)
    data$target4 = sub("-active", "", data$target4)
    plot_data8 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data8$time = 'TT8D'
  }
  
  # TT9D
  if(TRUE){
    data = read.table('../revised_data/TT9D/TT9D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = 'SC (c1)',
                           "2" = 'SC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-endosperm (c4)',
                           "5" = 'CZSC (c5)',
                           "6" = 'SC (c6)',
                           "7" = 'SC-ii (c7)',
                           "8" = 'SC-oi (c8)',
                           "9" = 'SC-oi (c9)',
                           "10" = 'SC-ii (c10)',
                           "11" = 'SC (c11)',
                           "12" = 'CZSC (c12)',
                           "13" = 'CZSC (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'CZSC (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'SC (c17)',
                           "18" = 'SC-ii (c18)',
                           "19" = 'SC-dividing-cell (c19)',
                           "20" = 'Endosperm-PCD (c20)',
                           "21" = 'CZSC (c21)',
                           "22" = 'CZSC-phloem-xylem (c22)',
                           "23" = 'SC (c23)',
                           "24" = 'EP-dividing-cell (c24)',
                           "25" = 'SC (c25)',
                           "26" = 'SC-SUS (c26)',
                           "27" = 'SC (c27)',
                           "28" = 'Endosperm-active (c28)')
    data$target2 <- recode(data$target, 
                           "1" = 'SC (c1)',
                           "2" = 'SC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC-endosperm (c4)',
                           "5" = 'CZSC (c5)',
                           "6" = 'SC (c6)',
                           "7" = 'SC-ii (c7)',
                           "8" = 'SC-oi (c8)',
                           "9" = 'SC-oi (c9)',
                           "10" = 'SC-ii (c10)',
                           "11" = 'SC (c11)',
                           "12" = 'CZSC (c12)',
                           "13" = 'CZSC (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'CZSC (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'SC (c17)',
                           "18" = 'SC-ii (c18)',
                           "19" = 'SC-dividing-cell (c19)',
                           "20" = 'Endosperm-PCD (c20)',
                           "21" = 'CZSC (c21)',
                           "22" = 'CZSC-phloem-xylem (c22)',
                           "23" = 'SC (c23)',
                           "24" = 'EP-dividing-cell (c24)',
                           "25" = 'SC (c25)',
                           "26" = 'SC-SUS (c26)',
                           "27" = 'SC (c27)',
                           "28" = 'Endosperm-active (c28)')
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = sub("-PCD", "", data$source3)
    data$target4 = sub("-PCD", "", data$target3)
    data$source4 = sub("-active", "", data$source4)
    data$target4 = sub("-active", "", data$target4)
    plot_data9 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data9$time = 'TT9D'
  }
  
  # TT11D
  if(TRUE){
    data = read.table('../revised_data/TT11D/TT11D_net_lr.txt', header = TRUE, sep = ',')
    data$source2 <- recode(data$source, 
                           "1" = 'CZSC (c1)',
                           "2" = 'SC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC (c4)',
                           "5" = 'SC-oi (c5)',
                           "6" = 'SC (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'CZSC (c9)',
                           "10" = 'SC (c10)',
                           "11" = 'CZSC (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'CZSC-phloem-xylem (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'SC (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'Endosperm (c17)')
    data$target2 <- recode(data$target, 
                           "1" = 'CZSC (c1)',
                           "2" = 'SC (c2)',
                           "3" = 'SC (c3)',
                           "4" = 'SC (c4)',
                           "5" = 'SC-oi (c5)',
                           "6" = 'SC (c6)',
                           "7" = 'SC (c7)',
                           "8" = 'SC (c8)',
                           "9" = 'CZSC (c9)',
                           "10" = 'SC (c10)',
                           "11" = 'CZSC (c11)',
                           "12" = 'SC (c12)',
                           "13" = 'CZSC-phloem-xylem (c13)',
                           "14" = 'SC (c14)',
                           "15" = 'SC (c15)',
                           "16" = 'SC (c16)',
                           "17" = 'Endosperm (c17)')
    data$source3 = sub(" .*", "", data$source2)
    data$target3 = sub(" .*", "", data$target2)
    data$source4 = sub("-PCD", "", data$source3)
    data$target4 = sub("-PCD", "", data$target3)
    data$source4 = sub("-active", "", data$source4)
    data$target4 = sub("-active", "", data$target4)
    plot_data10 = data[data$interaction_name_2 %in% c('ALE1->GSO2','RALFL33->AGP16','RALFL33->AGP20','LP1->CAM1','CYSB->MSS1'),]
    plot_data10$time = 'TT11D'
  }
  
  plot_data = rbind(plot_data1, plot_data2, plot_data3, plot_data4, plot_data5,
                    plot_data6, plot_data7, plot_data8, plot_data9, plot_data10)
  plot_data$treat = substr(plot_data$time,1,2)
  plot_data$sourceTarget = paste(plot_data$source4, '->', plot_data$target4, sep = '')
  plot_data$time = factor(plot_data$time, levels = c('CK3D', 'CK5D', 'CK7D', 'CK9D', 'CK11D','TT5D','TT7D','TT8D','TT9D','TT11D'))
  
  plot_data$interaction_name_3 = plot_data$interaction_name_2
  plot_data$interaction_name_3 = gsub('AT1G47710','Bch10G007420',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AT3G16920','Bch08G039110',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AT3G52370','Bch04G008210',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AtCDC48B','CDC48B',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AtCDC48C','CDC48C',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AT1G66250','Bch02G019130',plot_data$interaction_name_3)
  plot_data$interaction_name_3 = gsub('AT3G13560','Bch05G039790',plot_data$interaction_name_3)
  
  ggplot(plot_data, aes(x = time, y = sourceTarget, size = prob, color =  pval)) +
    geom_point(alpha = 0.6) + 
    scale_size_continuous(range = c(2,10)) + 
    scale_color_gradient(low = '#228B22', high = '#DAA520') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    # coord_flip()+
    facet_grid(~plot_data$interaction_name_3, scales = 'free', space = 'free') +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          text = element_text(size = 15),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'white'),
          plot.background = element_rect(fill="white"))
  
  ggsave('../revised_data/cellchat_pathway.svg', plot = last_plot(), height = 16, width = 12)
  ggsave('../CK_batch3_10samples/CK_cellchat_pathway2_v2.png', plot = last_plot(), height = 10, width = 9, bg = "transparent")
  
  
}
