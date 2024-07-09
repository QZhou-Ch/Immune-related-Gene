rm(list = ls())
library(doParallel)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(org.Hs.eg.db)
library(GEOquery)
library(sva)
library(affy)
library(limma)
library(DESeq2)
library(glmnet)
library(TissueEnrich)
library(clusterProfiler)
library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(survivalROC)
library(factoextra)
library(maftools)
library(RColorBrewer)
library(gghalves)
library(ggsci)
library(ggalluvial)
library(ggpubr)
library(ggridges)
library(ComplexHeatmap)
library(cowplot)
library(magrittr)
library(furrr)





library(tidyverse)

#####Figure2#####
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2A_D.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2E.rdata')


F2A_GPL570_KM_Train_Plot <- ggsurvplot(GPL570_KMfit_Train, data = GPL570_Train_clust, 
                                       risk.table = T, 
                                       conf.int = T,
                                       pval = T, pval.size = 6,
                                       #surv.median.line = "hv",#添加中位生存线
                                       xlab = "Follow up time(DAY)",font.x = 25,
                                       ylab = "OS Probability",font.y = 25,
                                       #legend = c(0.8,0.75),# 指定图例位置
                                       #font.tickslab = c(18, "plain", "darkgreen"),
                                       xlim = c(0,4000),break.x.by = 1000)
F2A_GPL570_KM_Train_Plot$plot <- F2A_GPL570_KM_Train_Plot$plot +
  xlab(NULL)+
  font("legend.title", size = 18)+ 
  font("legend.text", size = 18)
F2A_GPL570_KM_Train_Plot$table <- F2A_GPL570_KM_Train_Plot$table + 
  scale_y_discrete(labels = c('High','Low'))+
  theme(axis.title.x = element_text(size =25),axis.title.y = element_text(size =25))
F2A_GPL570_KM_Train_Plot <- plot_grid(F2A_GPL570_KM_Train_Plot$plot,F2A_GPL570_KM_Train_Plot$table,
                                      nrow = 2,rel_heights = c(3,1))


F2B_GPL570_Train_ROC_Plot <- ggplot(GPL570_Train_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
                                         round(GPL570_Train_clust_ROC$AUC
                                               [which(GPL570_Train_clust_ROC$Time_point == time_name[1])][1],3)),
                                   paste("AUC of 5 year =",
                                         round(GPL570_Train_clust_ROC$AUC
                                               [which(GPL570_Train_clust_ROC$Time_point == time_name[2])][1],3)),
                                   paste("AUC of 10 year =",
                                         round(GPL570_Train_clust_ROC$AUC
                                               [which(GPL570_Train_clust_ROC$Time_point == time_name[3])][1],3)))) +
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))

plot_F2A_B <- plot_grid(F2A_GPL570_KM_Train_Plot,F2B_GPL570_Train_ROC_Plot,
                        ncol = 2,align = "v",
                        labels = c("A", "B"))

a = 3
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F2_AB.pdf',width=a*6, height=a*3)
plot_F2A_B
dev.off()



F2C_GPL570_KM_Test_Plot <- ggsurvplot(GPL570_KMfit_Test, data = GPL570_Test_clust, 
                                      risk.table = T, 
                                      conf.int = T,
                                      pval = T, pval.size = 6,
                                      #surv.median.line = "hv",#添加中位生存线
                                      xlab = "Follow up time(DAY)",font.x = 25,
                                      ylab = "OS Probability",font.y = 25,
                                      #legend = c(0.8,0.75),# 指定图例位置
                                      #font.tickslab = c(18, "plain", "darkgreen"),
                                      xlim = c(0,4000),break.x.by = 1000)
F2C_GPL570_KM_Test_Plot$plot <- F2C_GPL570_KM_Test_Plot$plot +
  xlab(NULL)+
  font("legend.title", size = 18)+ 
  font("legend.text", size = 18)
F2C_GPL570_KM_Test_Plot$table <- F2C_GPL570_KM_Test_Plot$table + 
  scale_y_discrete(labels = c('High','Low'))+
  theme(axis.title.x = element_text(size =25),axis.title.y = element_text(size =25))
F2C_GPL570_KM_Test_Plot <- plot_grid(F2C_GPL570_KM_Test_Plot$plot,F2C_GPL570_KM_Test_Plot$table,
                                     nrow = 2,rel_heights = c(3,1))



F2D_GPL570_Test_ROC_Plot <- ggplot(GPL570_Test_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
                                         round(GPL570_Test_clust_ROC$AUC
                                               [which(GPL570_Test_clust_ROC$Time_point == time_name[1])][1],3)),
                                   paste("AUC of 5 year =",
                                         round(GPL570_Test_clust_ROC$AUC
                                               [which(GPL570_Test_clust_ROC$Time_point == time_name[2])][1],3)),
                                   paste("AUC of 10 year =",
                                         round(GPL570_Test_clust_ROC$AUC
                                               [which(GPL570_Test_clust_ROC$Time_point == time_name[3])][1],3)))) +
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))


plot_F2C_D <- plot_grid(F2C_GPL570_KM_Test_Plot,F2D_GPL570_Test_ROC_Plot, 
                        ncol = 2,align = "v",
                        labels = c("C","D"))
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F2_CD.pdf',width=a*6, height=a*3)
plot_F2C_D
dev.off()

F2E_GPL570_estimate_Plot <- ggplot(GPL570_estimate, aes(x = Group, y = value,fill = Group))+
  geom_boxplot()+
  stat_compare_means(method = 'wilcox.test',size = 4.2, vjust = -0.5,hjust = -0.5)+
  facet_grid(.~ESTIMATE)+
  labs(y = "ESTIMATE", x = "")+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        axis.text.x = element_text(size = 15),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        strip.text = element_text(size = 15))
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F2_E.pdf',width=a*6, height=a*3)
F2E_GPL570_estimate_Plot
dev.off()


#####Figure3#####

load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F3A_F.rdata')


F3A_TCGA_KM_OS_Plot <- ggsurvplot(TCGA_fit_OS, data = tcga_clust, 
                                  risk.table = T, 
                                  conf.int = T,
                                  pval = T, pval.size = 6,
                                  #surv.median.line = "hv",#添加中位生存线
                                  xlab = "Follow up time(DAY)",font.x = 25,
                                  ylab = "OS Probability",font.y = 25,
                                  #legend = c(0.8,0.75),# 指定图例位置
                                  #font.tickslab = c(18, "plain", "darkgreen"),
                                  xlim = c(0,4000),break.x.by = 1000)

F3A_TCGA_KM_OS_Plot$plot <- F3A_TCGA_KM_OS_Plot$plot +
  xlab(NULL)+
  font("legend.title", size = 18)+ 
  font("legend.text", size = 18)
F3A_TCGA_KM_OS_Plot$table <- F3A_TCGA_KM_OS_Plot$table + 
  scale_y_discrete(labels = c('High','Low'))+
  theme(axis.title.x = element_text(size =25),axis.title.y = element_text(size =25))
F3A_TCGA_KM_OS_Plot <- plot_grid(F3A_TCGA_KM_OS_Plot$plot,F3A_TCGA_KM_OS_Plot$table,
                                      nrow = 2,rel_heights = c(3,1))


F3B_TCGA_OS_ROC_Plot <- ggplot(TCGA_OS_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
                                         round(TCGA_OS_clust_ROC$AUC
                                               [which(TCGA_OS_clust_ROC$Time_point == time_name[1])][1],3)),
                                   paste("AUC of 5 year =",
                                         round(TCGA_OS_clust_ROC$AUC
                                               [which(TCGA_OS_clust_ROC$Time_point == time_name[2])][1],3)),
                                   paste("AUC of 10 year =",
                                         round(TCGA_OS_clust_ROC$AUC
                                               [which(TCGA_OS_clust_ROC$Time_point == time_name[3])][1],3)))) +
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))


F3C_TCGA_KM_PFI_Plot <- ggsurvplot(TCGA_fit_PFI, data = PFI_clust, 
                                   risk.table = T, 
                                   conf.int = T,
                                   pval = T, pval.size = 6,
                                   #surv.median.line = "hv",#添加中位生存线
                                   xlab = "Follow up time(DAY)",font.x = 25,
                                   ylab = "PFS Probability",font.y = 25,
                                   #legend = c(0.8,0.75),# 指定图例位置
                                   #font.tickslab = c(18, "plain", "darkgreen"),
                                   xlim = c(0,4000),break.x.by = 1000)

F3C_TCGA_KM_PFI_Plot$plot <- F3C_TCGA_KM_PFI_Plot$plot +
  xlab(NULL)+
  font("legend.title", size = 18)+ 
  font("legend.text", size = 18)
F3C_TCGA_KM_PFI_Plot$table <- F3C_TCGA_KM_PFI_Plot$table + 
  scale_y_discrete(labels = c('High','Low'))+
  theme(axis.title.x = element_text(size =25),axis.title.y = element_text(size =25))
F3C_TCGA_KM_PFI_Plot <- plot_grid(F3C_TCGA_KM_PFI_Plot$plot,F3C_TCGA_KM_PFI_Plot$table,
                                 nrow = 2,rel_heights = c(3,1))

F3D_TCGA_PFI_ROC_Plot <- ggplot(TCGA_PFI_clust_ROC,aes(x = False_positive,y = True_positive,colour = Time_point)) +
  geom_line(size = 0.5) +
  guides(colour = guide_legend(title = NULL))+
  scale_colour_discrete(labels = c(paste("AUC of 3 year =",
                                         round(TCGA_PFI_clust_ROC$AUC
                                               [which(TCGA_PFI_clust_ROC$Time_point == time_name[1])][1],3)),
                                   paste("AUC of 5 year =",
                                         round(TCGA_PFI_clust_ROC$AUC
                                               [which(TCGA_PFI_clust_ROC$Time_point == time_name[2])][1],3)),
                                   paste("AUC of 10 year =",
                                         round(TCGA_PFI_clust_ROC$AUC
                                               [which(TCGA_PFI_clust_ROC$Time_point == time_name[3])][1],3)))) +
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))


F3E_TCGA_OS_forest_Plot <- ggforest(TCGA_OS_model,  #coxph得到的Cox回归结果
                                    data = TCGA_multiCOX_PFI,  #数据集
                                    #标题
                                    main = NULL,  
                                    #前三列距离
                                    cpositions = c(0.05, 0.15, 0.35),  
                                    #字体大小
                                    fontsize = 1, 
                                    #相对变量的数值标签，也可改为1
                                    refLabel = 'reference', 
                                    #保留HR值以及95%CI的小数位数
                                    noDigits = 3)+
  theme_void()+
  ggtitle('Hazard ratio of BRCA (OS)')+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        plot.title = element_text(size = 25,hjust = 0.5),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))



F3F_TCGA_PFI_forest_Plot <- ggforest(TCGA_PFI_model,  #coxph得到的Cox回归结果
                                     data = TCGA_multiCOX_PFI,  #数据集
                                     #标题
                                     main = NULL,  
                                     #前三列距离
                                     cpositions = c(0.05, 0.15, 0.35),  
                                     #字体大小
                                     fontsize = 1, 
                                     #相对变量的数值标签，也可改为1
                                     refLabel = 'reference', 
                                     #保留HR值以及95%CI的小数位数
                                     noDigits = 3)+
  theme_void()+
  ggtitle('Hazard ratio of BRCA (PFS)')+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        plot.title = element_text(size = 25,hjust = 0.5),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))


a = 3
plot_F3A_B <- plot_grid(F3A_TCGA_KM_OS_Plot, F3B_TCGA_OS_ROC_Plot, 
                        ncol = 2,align = "v",
                        labels = c("A", "B"))
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F3_AB.pdf',width=a*6, height=a*3)
plot_F3A_B
dev.off()

plot_F3C_D <- plot_grid(F3C_TCGA_KM_PFI_Plot,F3D_TCGA_PFI_ROC_Plot,
                        ncol = 2,align = "v",
                        labels = c("C","D"))
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F3_CD.pdf',width=a*6, height=a*3)
plot_F3C_D
dev.off()

plot_F3E_F <- plot_grid(F3E_TCGA_OS_forest_Plot,F3F_TCGA_PFI_forest_Plot,
                        ncol = 2,align = "v",
                        labels = c("E","F"))

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F3_EF.pdf',width=a*6, height=a*3)
plot_F3E_F
dev.off()


#####Figure4#####
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F4A_D.rdata')


F4A_Immune_Subtype_Box_plot <- ggplot(Immune_Subtype , aes(x = Immune.Subtype, y = Immune.Score, fill = Immune.Subtype))+
  geom_boxplot(aes(x = Immune.Subtype,y = Immune.Score, fill = Immune.Subtype),
               outlier.shape = NA,
               width = .07,
               color = "black")+
  stat_compare_means(label.y = 3.5,size = 5)+
  stat_compare_means(comparisons = my_comparisons)+
  geom_half_violin(aes(fill = Immune.Subtype),
                   position = position_nudge(x = .15, y = 0),
                   width = .2,
                   adjust=1.5, trim=FALSE, colour=NA, side = 'r') +
  geom_point(aes(x = as.numeric(Immune.Subtype)-0.1,
                 y = Immune.Score,color = Immune.Subtype),
             position = position_jitter(width = .05),size = .25, shape = 20) +
  scale_color_jco() +
  scale_fill_jco() +
  theme_minimal()+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        axis.text.x = element_text(size = 15),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))

a = 3
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_A.pdf',width=16,height=5)
F4A_Immune_Subtype_Box_plot
dev.off()


F4B_RNASeq_FPKM_Plot <- ggplot(RNASeq_FPKM, aes(x = Group, y = value,fill = Group))+
  geom_boxplot()+
  stat_compare_means(method = 'wilcox.test',size = 4, vjust = 2,hjust = 0.5)+
  facet_grid(.~FPKM)+
  labs(y = "FPKM", x = "",size = 15)+
  #scale_x_discrete(limits = c(`PD-1`,`PD-L1`,`CTLA4`))+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        axis.title.x = element_text(size =25),
        axis.title.y = element_text(size =25),
        axis.text.x = element_text(size = 15),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        strip.text = element_text(size = 12))

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_B.pdf',width=10,height=5)
F4B_RNASeq_FPKM_Plot
dev.off()



F4C1_Immune_Sankey_plot <- ggplot(as.data.frame(Sankey_Immune),
                                  aes(axis1 = Immune_Subtype, axis2 = Group,y = count)) +
  scale_x_discrete(limits = c("Immune_Subtype", "Group"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Group),alpha = 0.8) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()+
  ylab("SubType") +
  theme(legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_C1.pdf',width=5,height=5)
F4C1_Immune_Sankey_plot
dev.off()


F4C2_IHC_Sankey_plot <- ggplot(as.data.frame(Sankey_IHC),
                               aes(axis1 = Group, axis2 = IHC_Subtype, y = count)) +
  scale_x_discrete(limits = c("Group", "IHC_Subtype"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Group),alpha = 0.8) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()+
  ylab("") +
  theme(legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_C2.pdf',width=5,height=5)
F4C2_IHC_Sankey_plot
dev.off()

F4D_Lollipop_plot <- ggdotchart(wilcox_result, x = "name", y = "dis",
                                color = 'p-value', palette = c("#8e9eab","#ffba08","#e85d04","#d00000"), # 修改颜色
                                size = 16, legend.title = 'Wilcox p-value',
                                label = 'H_L', font.label = list(color = "white", size = 16, 
                                                                 vjust = 0.5), 
                                add = "segments", add.params = list(color = "lightgray", size = 3),
                                #ggtheme = theme_pubr(),
                                xlab = "", ylab = '',
                                sorting = "descending", group = "group", 
                                rotate = TRUE) +
  scale_y_continuous(breaks = c(-0.05,0,0.05),labels = c('Low_Score','|','High_Score'))+
  scale_x_discrete(position = 'top')+
  theme(legend.position = c(0.90,0.08), 
        #axis.text=element_text(size=18), 
        axis.ticks.x = element_blank())

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_D.pdf',width=10,height=20)
F4D_Lollipop_plot
dev.off()




#F4E
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F4_E.pdf',width=16,height=5)
oncoplot(maf = brca_maf,top = 20,
         clinicalFeatures = c('Group','IHC.Subtype','Immune.Subtype'), sortByAnnotation = T,
         annotationColor = list(Group = c('High_Score'="#d03027",'Low_Score'="#004977"),
                                IHC.Subtype = c('HR-/HER2+' = '#7ac70c', 'HR+/HER2-' = '#faa918',
                                                'HR+/HER2+'='#1cb0f6','NA' = '#d0d2d3','TNBC' = '#8549ba' ),
                                Immune.Subtype = c('C1'='#8ee000','C2'='#ffc715',
                                                   'C3'='#14d4f4','C4'= '#cfcfcf','C6'= '#a560e8')),
         leftBarData = focal_CNV_per,leftBarLims = c(0,20),gene_mar = 3,
         fontSize = 0.5,logColBar = T,
         legend_height = 3,
         legendFontSize = 1, annotationFontSize = 1, titleFontSize = 1.5)
dev.off()



#####Figure5#####
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F5A_E.rdata')


F5A_TCGA_DE_DOT_Plot <- ggscatter(TCGA_DE_DOT,x = "log2FoldChange",y = 'logp',#指定数据集，x轴，y轴
                                  color = 'Group',#颜色按照差异低表达，无差异表达，差异高表达分组
                                  palette = c("#619cff","#d7d7d8","#f8766d"),#指定分组对应的颜色
                                  size = 3, gamma = 0.6,
                                  #shape = 'Type',#形状按照mRNA，lncRNA，miRNA
                                  #alpha = "Expression",#透明度按照表达量展示
                                  #size=DotPlot$Expression,#指定点的大小
                                  #label = DotPlot$genename,
                                  #label.select = DotPlot$genename[!is.na(DotPlot$genename)],#展示前10个差异基因的基因名
                                  #font.label = 20,
                                  legend.title = NULL,
                                  repel = T)+
  theme_classic2()+
  geom_hline(yintercept=1.30,linetype ="dashed")+
  geom_vline(xintercept=c(-1,1),linetype ="dashed")+
  labs(title="Different Expression Between \nHigh-Score and Low-Score Group",x="Log2 FoldChange", y="-Log10(Adj P-value)")+
  font("title", size = 30)+
  #guides(fill = guide_legend(title = NULL))+
  #font("legend.text", size = 30)+
  # #增加一个图层
  # ggrepel::geom_label_repel(aes(label = genename),data = TCGA_DE_DOT[-which(is.na(TCGA_DE_DOT$genename)),],
  #                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100),)+
  theme(legend.title = element_text(size = 0),
        legend.text = element_text(size = 30),
        legend.key.size = unit(45, "pt"),
    plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size =30),
        axis.title.y = element_text(size = 30))

F5A_TCGA_DE_DOT_Plot

a = 4
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F5_A.pdf',width=a*3, height=a*3)
F5A_TCGA_DE_DOT_Plot
dev.off()



F5B_immune_cell_GSEA_Plot <- enrichplot::gseaplot2(immune_cell,1:4,ES_geom= "line",base_size = 25)
F5B_immune_cell_GSEA_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F5_B.pdf',width=a*3, height=a*3)
F5B_immune_cell_GSEA_Plot
dev.off()



F5C_enrichGO_Plot <- dotplot(ego, color = "p.adjust", title = "GO Enrichment",
                        orderBy = "x",showCategory= 20,font.size	= 12)+ 
  facet_grid(ONTOLOGY~.,scale="free",space="free")+
  font("xlab", size = 20)+
  font("title", size = 25)

F5D_enrichKEGG_Plot <- dotplot(ekegg,color = "p.adjust", title = "KEGG Enrichment", 
                               orderBy = "x",showCategory= 20,font.size	= 12)+
  font("xlab", size = 20)+
  font("title", size = 25)

plot_F5C_D <- plot_grid(F5C_enrichGO_Plot,F5D_enrichKEGG_Plot,
                        ncol = 2,rel_widths = c(0.5,0.5),
                        labels = c("C","D"))
a = 3
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F5_CD.pdf',width=a*7, height=a*3)
plot_F5C_D
dev.off()

F5E_heatmap_Plot <- Heatmap(MSigDB_Heat,
                            heatmap_width = unit(14, "cm"), heatmap_height = unit(24, "cm"),
                            #width = unit(300, "pt"), height = unit(750, "pt"),
                            row_gap = unit(2, "mm"),
                            heatmap_legend_param = list(title = "NES"),
                            col = circlize::colorRamp2(c(-3, 0, 3), c("#004977", "white", "#d03027")),
                            cluster_rows = T, cluster_columns = F,
                            row_names_gp = gpar(fontsize = 10),
                            show_row_dend = F,
                            row_km = 2, row_title = NULL,
                            show_column_names = F, row_names_side = "left",row_dend_side = "right",
                            bottom_annotation = HeatmapAnnotation(Score = as.numeric(ANNO_Heat$score),
                                                                  col = list(Score = circlize::colorRamp2(c(min(as.numeric(ANNO_Heat$score)),
                                                                                                            max(as.numeric(ANNO_Heat$score))),
                                                                                                          c('#ffffff','#338d11')))),
                            top_annotation = HeatmapAnnotation(Group = factor(ANNO_Heat$Group,levels = c('High_Score','Low_Score')),
                                                               col = list(Group = c('High_Score'="#d03027",'Low_Score'="#004977"))))

F5E_heatmap_Plot

a = 9
pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/ToAI/F5_E.pdf',width=a*1, height=a*2)
F5E_heatmap_Plot
dev.off()

























