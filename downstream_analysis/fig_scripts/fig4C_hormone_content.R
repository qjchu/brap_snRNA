library(ggplot2)
library(dplyr)

getwd()
setwd('scripts/')

#################### IAA ng/mL
if(TRUE){
  data <- data.frame(
    Group = rep(c("CK6D", "CK10D", "TT6D", "TT10D"), each = 3),
    Value = c(79.82,79.13,78.25, 151.47,151.64,150.86,
              34.29,34.34,34.19, 26.64,7.27,7.43)
  )
  
  summary_data <- data %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value)
    )
  
  print(summary_data)
  
  t.test(c(79.82,79.13,78.25), c(34.29,34.34,34.19))
  t.test(c(151.47,151.64,150.86), c(26.64,7.27,7.43))
  
  summary_data$Group = factor(summary_data$Group, levels = c('CK6D','TT6D','CK10D','TT10D'))
  ggplot(summary_data, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.7) +  
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD), 
      width = 0.2, 
      color = "black",
      size = 0.8
    ) +
    scale_fill_manual(values = rep(c( '#556B2F', '#DAA520'),3)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill="transparent"),
          legend.position = "none") 
  
  ggsave('../revised_data/CK_TT_geneset/CK_TT_IAA_content_v3.svg', plot = last_plot(), height = 3.5, width = 5, bg = "transparent")

}


#################### GA3 ng/mL
if(TRUE){
  data <- data.frame(
    Group = rep(c("CK6D", "CK10D", "TT6D", "TT10D"), each = 3),
    Value = c(0.14, 0.11, 0.13, 0.049, 0.05, 0.051,
              0.06, 0.06, 0.07, 0.07, 0.01, 0.03 )
  )
  
  summary_data <- data %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value)
    )
  
  print(summary_data)
  
  t.test(c(0.14, 0.11, 0.13), c(0.06, 0.06, 0.07))
  t.test(c(0.049, 0.05, 0.051), c(0.07, 0.01, 0.03))
  
  summary_data$Group = factor(summary_data$Group, levels = c('CK6D','TT6D','CK10D','TT10D'))
  ggplot(summary_data, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.7) +  
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD),  
      width = 0.2, 
      color = "black",
      size = 0.8
    ) + 
    scale_fill_manual(values = rep(c('#556B2F',  '#DAA520'),3)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill="transparent"),
          legend.position = "none") 
  
  ggsave('../revised_data/CK_TT_geneset/CK_TT_GA3_content_v3.svg', plot = last_plot(), height = 3.5, width = 5, bg = "transparent")
  
}

#################### ABA ng/mL
if(TRUE){
  data <- data.frame(
    Group = rep(c("CK6D", "CK10D", "TT6D", "TT10D"), each = 3),
    Value = c(74.32,74.35,74.16, 31.02,30.68,31.22,
              148.51,148.62,151.41,364.42,319.60,480.50)
  )
  
  summary_data <- data %>%
    group_by(Group) %>%
    summarise(
      Mean = mean(Value),
      SD = sd(Value)
    )
  
  print(summary_data)
  
  t.test(c(74.32,74.35,74.16), c(148.51,148.62,151.41))
  t.test(c(31.02,30.68,31.22), c(364.42,319.60,480.50))
  
  summary_data$Group = factor(summary_data$Group, levels = c('CK6D','TT6D','CK10D','TT10D'))
  ggplot(summary_data, aes(x = Group, y = Mean, fill = Group)) +
    geom_col(width = 0.7) + 
    geom_errorbar(
      aes(ymin = Mean - SD, ymax = Mean + SD), 
      width = 0.2, 
      color = "black",
      size = 0.8
    ) + 
    scale_fill_manual(values = rep(c(  '#556B2F','#DAA520'),3)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 8),
          axis.title = element_blank(),
          text = element_text(size = 12),
          axis.line = element_line(color = 'darkgray'),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill="transparent"),
          legend.position = "none") 
  
  ggsave('../revised_data/CK_TT_geneset/CK_TT_ABA_content_v3.svg', plot = last_plot(), height = 3.5, width = 5, bg = "transparent")
  

}
