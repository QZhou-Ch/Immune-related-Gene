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
library(magrittr)
library(furrr)





library(tidyverse)


load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/4.1TCGA_Score.rdata')



#####5个免疫亚型的评分，进行kruskal检验，并绘制箱式图########################################

Immune_Subtype <- tcga_pscore[,c(3,13)] %>% 
  filter(!is.na(Immune.Subtype)) %>% 
  mutate(Immune.Subtype = factor(Immune.Subtype),
         Immune.Score = as.numeric(score)) %>% 
  select(Immune.Subtype,Immune.Score)

#云雨图不好看


my_comparisons <- list( c("C1", "C2"), c("C1", "C3"), c("C1", "C4"), c("C1", "C6"),
                        c("C2", "C3"), c("C2", "C4"), c("C2", "C6"),
                        c("C3", "C4"), c("C3", "C6"),
                        c("C4", "C6"))

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
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))

F4A_Immune_Subtype_Box_plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_A.pdf',width=15,height=15)
F4A_Immune_Subtype_Box_plot
dev.off()



#####3个RNASeq-FPKM在两个分组中的boxplot#######################################################
##PD-L1 CD274 ENSG00000120217
##PD-1 CD279 ENSG00000188389
##CTLA4 CD152 ENSG00000163599
##LAG3 CD223  ENSG00000089692
##TIM-3 CD366 ENSG00000135077
RNASeq_FPKM <- tcga_rnaseq_FPKM %>% 
  as.data.frame() %>% 
  dplyr::select(c('sample','ENSG00000120217','ENSG00000188389',
                  'ENSG00000163599')) %>% 
  left_join(tcga_pscore[,c(2:4)], by = 'sample') %>% 
  filter(!is.na(score)) %>% 
  as_tibble() %>% 
  dplyr::select('sample','Group',
                `PD-L1` = 'ENSG00000120217',
                `PD-1` = 'ENSG00000188389',
                CTLA4 = 'ENSG00000163599') %>% 
  #转换成长数据
  pivot_longer(-c(sample,Group),names_to = "FPKM", values_to = "value") %>% 
  mutate(value = as.numeric(value))

RNASeq_result <- tcga_rnaseq_FPKM %>% 
  as.data.frame() %>% 
  dplyr::select(c('sample','ENSG00000120217','ENSG00000188389',
                  'ENSG00000163599')) %>% 
  left_join(tcga_pscore[,c(2:4)], by = 'sample') %>% 
  filter(!is.na(score)) %>% 
  as_tibble() %>% 
  dplyr::select('sample','score',
                `PD-L1` = 'ENSG00000120217',
                `PD-1` = 'ENSG00000188389',
                CTLA4 = 'ENSG00000163599')

cor_RNASeq_result <- map(RNASeq_result[,-c(1,2)],function(x){
  #x = immune_factor[,-c(1,2)][,18]
  cor_res = cor.test(x,as.numeric(RNASeq_result$score),method = "spearman",exact=FALSE)
  a = cor_res[["p.value"]]
  b = cor_res[["estimate"]][["rho"]]
  re = c(a,b)
  names(re) = c("pvalue","rho")
  return(re)
}) %>% 
  reduce(rbind) %>% 
  as.data.frame() 

rownames(cor_RNASeq_result) <- c('PD-L1', 'PD-1', 'CTLA4')


RNASeq_result <- ggplot(data = RNASeq_result, aes(x = 'score', y = 'CTLA4'))+
  geom_point(color="red")+
  stat_smooth(score~CTLA4, method="lm", se=T)+
  stat_cor(data = RNASeq_result, method = "spearman",exact=FALSE)


#stat_cor(data=dat, method = "pearson")意为用pearson相关进行相关性分析，可以自行更改方法





##做图

F4B_RNASeq_FPKM_Plot <- ggplot(RNASeq_FPKM, aes(x = Group, y = value,fill = Group))+
  geom_boxplot()+
  stat_compare_means(method = 'wilcox.test',size = 4, vjust = 2,hjust = 0.5)+
  facet_grid(.~FPKM)+
  labs(y = "FPKM", x = "",size = 15)+
  #scale_x_discrete(limits = c(`PD-1`,`PD-L1`,`CTLA4`))+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        strip.text = element_text(size = 12))
F4B_RNASeq_FPKM_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_B.pdf',width=15,height=15)
F4B_RNASeq_FPKM_Plot
dev.off()



#####高分低分组，与PAM50和IHC分类的关系，并绘Sankey图##########################################

Sankey <- tcga_pscore[,c(4,12:13)] %>% 
  #filter(!((is.na(IHC.Subtype)) | (IHC.Subtype == 'NA') | (is.na(Immune.Subtype)))) %>% 
  rename(IHC_Subtype = IHC.Subtype,
         Immune_Subtype = Immune.Subtype) %>% 
  select(Group,IHC_Subtype,Immune_Subtype)




#####做图1
library(ggforce)
Sankey_plot <- Sankey %>% 
  filter(!(IHC_Subtype=='NA'| Immune_Subtype=="NA")) %>% 
  mutate(IHC = factor(IHC_Subtype),Group = factor(Group),Immune = factor(Immune_Subtype)) %>% 
  select('IHC','Group','Immune','count') %>% 
  as.data.frame() %>% gather_set_data(1:3) %>% 
  ggplot(aes(x, id = id, split = y, value = count)) +
  geom_parallel_sets(aes(fill = Group), alpha = 0.3, axis.width = 0.1) +
  scale_x_discrete(limits = c('IHC','Group','Immune'))+
  geom_parallel_sets_axes(axis.width = 0.1) +
  geom_parallel_sets_labels(colour = 'white')
Sankey_plot




#####做图2
Sankey_IHC <- Sankey %>% 
  select(Group,IHC_Subtype) %>% 
  filter(!(IHC_Subtype == 'NA')) %>% 
  group_by(Group,IHC_Subtype) %>% 
  summarise(., count = n()) %>% 
  as_tibble()

IHC_hspart <- data.frame(IHC_Subtype = c('HR-/HER2+','HR+/HER2-','HR+/HER2+','TNBC'),
                         partition = c((Sankey_IHC$count[1]/(Sankey_IHC$count[1]+Sankey_IHC$count[5])),
                                       (Sankey_IHC$count[2]/(Sankey_IHC$count[2]+Sankey_IHC$count[6])),
                                       (Sankey_IHC$count[3]/(Sankey_IHC$count[3]+Sankey_IHC$count[7])),
                                       (Sankey_IHC$count[4]/(Sankey_IHC$count[4]+Sankey_IHC$count[8]))))

IHC_hspart

Sankey_Immune <- Sankey %>% 
  select(Group,Immune_Subtype) %>% 
  filter(!(is.na(Immune_Subtype))) %>% 
  group_by(Group,Immune_Subtype) %>% 
  summarise(., count = n()) %>% 
  as_tibble()

Immune_hspart <- data.frame(IHC_Subtype = c('C1','C2','C3','C4','C6'),
                            partition = c((Sankey_Immune$count[1]/(Sankey_Immune$count[1]+Sankey_Immune$count[6])),
                                       (Sankey_Immune$count[2]/(Sankey_Immune$count[2]+Sankey_Immune$count[7])),
                                       (Sankey_Immune$count[3]/(Sankey_Immune$count[3]+Sankey_Immune$count[8])),
                                       (Sankey_Immune$count[4]/(Sankey_Immune$count[4]+Sankey_Immune$count[9])),
                                       (Sankey_Immune$count[5]/(Sankey_Immune$count[5]+Sankey_Immune$count[10]))))
Immune_hspart


F4C1_Immune_Sankey_plot <- ggplot(as.data.frame(Sankey_Immune),
                             aes(axis1 = Immune_Subtype, axis2 = Group,y = count)) +
  scale_x_discrete(limits = c("Immune_Subtype", "Group"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Group)) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()+
  ylab("SubType") +
  theme(legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))
F4C1_Immune_Sankey_plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_C1.pdf',width=15,height=15)
F4C1_Immune_Sankey_plot
dev.off()


F4C2_IHC_Sankey_plot <- ggplot(as.data.frame(Sankey_IHC),
                          aes(axis1 = Group, axis2 = IHC_Subtype, y = count)) +
  scale_x_discrete(limits = c("Group", "IHC_Subtype"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Group)) +
  geom_stratum() + 
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_void()+
  ylab("") +
  theme(legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))
F4C2_IHC_Sankey_plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_C2.pdf',width=15,height=15)
F4C2_IHC_Sankey_plot
dev.off()

#####高分低分组，结合免疫相关分数，进行两组检验，并绘制箱线图################################

immune_factor <- tcga_pscore[,c(3:4,14:18,24:31)] %>% 
  filter(!is.na(Group)) %>% 
  select("score","Group",
#Tumor Immune Infiltrate=      
  "Leukocyte.Fraction","Stromal.Fraction","TIL.Regional.Fraction","Proliferation",
#Correlates of Somatic Variation = 
"Number.of.Segments","Aneuploidy.Score","Homologous.Recombination.Defects","Intratumor.Heterogeneity",
#Immunogenicity
"SNV.Neoantigens","Indel.Neoantigens",
#Genomic State
  "Silent.Mutation.Rate","Nonsilent.Mutation.Rate","Fraction.Altered",
#Adaptive Immune Receptor 
  #"BCR.Evenness","BCR.Shannon","BCR.Richness",
  #"TCR.Shannon","TCR.Richness","TCR.Evenness","CTA.Score"
  )

#####wilcox检验及相关分析

wilcox_result <- map(immune_factor[,-c(1,2)],function(x){
  #x = immune_factor[,-c(1,2)][,18]
  x = (x - min(x,na.rm = T))/(max(x,na.rm = T) - min(x,na.rm = T))
  wilcox_res = wilcox.test(x~immune_factor$Group, exact = F)
  a = wilcox_res[["p.value"]]
  b = wilcox_res[["statistic"]]
  c = cbind(x,immune_factor$Group) %>% 
    as.data.frame() %>% 
    group_by(V2) %>%
    summarize(median = median(as.numeric(x),na.rm = T))
  re = c(a,b,c$median)
  names(re) = c("pvalue","W","Median_High","Median_Low")
  return(re)
}) %>% 
  reduce(rbind) %>% 
  as.data.frame() %>% 
  mutate(res = ifelse((as.numeric(.$Median_High) - as.numeric(.$Median_Low))>0,'High','Low')) %>% 
  select("Median_High","Median_Low","W","pvalue","res")

rownames(wilcox_result) <- colnames(immune_factor[,c(3:ncol(immune_factor))])


cor_result <- map(immune_factor[,-c(1,2)],function(x){
  #x = immune_factor[,-c(1,2)][,18]
  cor_res = cor.test(x,as.numeric(immune_factor$score),method = "spearman", exact = F)
  a = cor_res[["p.value"]]
  b = cor_res[["estimate"]][["rho"]]
  re = c(a,b)
  names(re) = c("pvalue","rho")
  return(re)
}) %>% 
  reduce(rbind) %>% 
  as.data.frame() 

rownames(cor_result) <- colnames(immune_factor[,c(3:ncol(immune_factor))])

#####做图1

wilcox_result$name <- rownames(wilcox_result)



wilcox_result$dis <- (wilcox_result$Median_High - wilcox_result$Median_Low)
wilcox_result['SNV.Neoantigens','dis'] <- wilcox_result['SNV.Neoantigens','dis']+0.005
wilcox_result['Indel.Neoantigens','dis'] <- wilcox_result['SNV.Neoantigens','dis']+0.005
wilcox_result['Silent.Mutation.Rate','dis'] <- wilcox_result['SNV.Neoantigens','dis']+0.005
wilcox_result['Nonsilent.Mutation.Rate','dis'] <- wilcox_result['SNV.Neoantigens','dis']+0.005

wilcox_result$H_L <- ifelse(wilcox_result$dis > 0,"H","L")

wilcox_result$`p-value` <- NA
wilcox_result[which(wilcox_result$pvalue > 0.05),'p-value'] <- 'NA'
wilcox_result[which(wilcox_result$pvalue < 0.05 & wilcox_result$pvalue > 0.01),'p-value'] <- '*'
wilcox_result[which(wilcox_result$pvalue < 0.01 & wilcox_result$pvalue > 0.001),'p-value'] <- '**'
wilcox_result[which(wilcox_result$pvalue < 0.001),'p-value'] <- '***'
wilcox_result$`p-value` <- factor(wilcox_result$'p-value',levels = c('NA','*','**','***'))

wilcox_result$group <- NA
wilcox_result[which(wilcox_result$name %in% 
                      c("Leukocyte.Fraction",
                        "Stromal.Fraction",
                        "TIL.Regional.Fraction",
                        "Proliferation")),'group'] <- 'Tumor Immune Infiltrate'
wilcox_result[which(wilcox_result$name %in% 
                      c("Number.of.Segments",
                        "Aneuploidy.Score",
                        "Homologous.Recombination.Defects",
                        "Intratumor.Heterogeneity")),'group'] <- 'Somatic Variation'
wilcox_result[which(wilcox_result$name %in% 
                      c("SNV.Neoantigens","Indel.Neoantigens")),'group'] <- 'Immunogenicity'
wilcox_result[which(wilcox_result$name %in% 
                      c("Silent.Mutation.Rate","Nonsilent.Mutation.Rate","Fraction.Altered")),'group'] <- 'Genomic State'
wilcox_result$group <- factor(wilcox_result$group,levels = c('Immunogenicity','Genomic State',
                                                             'Somatic Variation','Tumor Immune Infiltrate'))


F4D_Lollipop_plot <- ggdotchart(wilcox_result, x = "name", y = "dis",
                             color = 'p-value', palette = c("#8e9eab","#ffba08","#e85d04","#d00000"), # 修改颜色
                             size = 6, legend.title = 'Wilcox p-value',
                             label = 'H_L', font.label = list(color = "white", size = 9, 
                                                              vjust = 0.5), 
                             add = "segments", add.params = list(color = "lightgray", size = 1.5),
                             #ggtheme = theme_pubr(),
                             xlab = "", ylab = '',
                             sorting = "descending", group = "group", 
                             rotate = TRUE) +
  scale_y_continuous(breaks = c(-0.05,0,0.05),labels = c('Low_Score','|','High_Score'))+
  theme(legend.position = c(0.90,0.18),  axis.ticks.x = element_blank())

F4D_Lollipop_plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_D.pdf',width=15,height=15)
F4D_Lollipop_plot
dev.off()



#####做图2
immune_factor_Box_plot <- immune_factor[,-1] %>% 
  pivot_longer(-c(Group),names_to = "Type", values_to = "value") %>% 
  split(.$Type) %>%
  map(~ggplot(., aes(x = Type, y = value,fill = Group)) + 
        geom_boxplot(outlier.alpha = 0.4,width = 0.2,outlier.size = 0.5)+
        stat_compare_means(method = 'wilcox.test', vjust = 0,hjust = -0.5,size = 4)+
        labs(y="",x="")+
        theme_minimal()+
        theme(panel.background=element_rect(colour="black",fill=NA),
              panel.grid=element_blank(),
              legend.position= c(0.9,0.1),
              legend.background=element_rect(colour=NA,fill=NA)))

patch <- stringr::str_c('F3_F_',names(immune_factor_Box_plot),".pdf")
pwalk(list(patch,immune_factor_Box_plot),ggsave,path = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/')






#####高分低分两组SNP的差异，并绘制瀑布图#######################################################

allmaf <- read.maf(maf = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.gz')

brca_clin <- as.data.frame(allmaf@variant.classification.summary) %>% 
  mutate(pat = str_sub(.[,1],start = 1,end = 12)) %>% 
  left_join(tcga_pscore[,c(1,4,11:13)], by = 'pat') %>% 
  filter(!is.na(Group)) %>% 
  #mutate(Group = factor(Group,ordered = T, levels = c('High_Score','Low_Score'))) %>% 
  #mutate(Group = factor(Group,levels = c('High_Score','Low_Score','High_Score'))) %>% 
  select(Tumor_Sample_Barcode,Group,PAM50.Subtype,IHC.Subtype,Immune.Subtype)


brca_maf <- read.maf(maf = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/TCGA.BRCA.varscan.6c93f518-1956-4435-9806-37185266d248.DR-10.0.somatic.maf.gz', 
                     clinicalData = brca_clin) %>% 
  subsetMaf(maf = ., tsb= brca_clin$Tumor_Sample_Barcode)

cnvs_per <- function(focal_CNV) {
  focal_CNV <- cbind(focal_CNV[,1], 
                     as.data.frame(apply(focal_CNV[,-1],c(1,2),function(x){
                       x = as.numeric(x)
                       #以0.3为cut-off值
                       if(x < -0.3){
                         x = 'DEL'
                       }else if(x > 0.3){
                         x = 'AMP'
                       }else{x = 'None'}
                       x = as.character(x)
                       return(x)
                     })))
  colnames(focal_CNV) <- str_replace_all(colnames(focal_CNV),pattern = "[.]",replacement = "-")
  rownames(focal_CNV) <- focal_CNV[,1]
  
  focal_CNV2 <- focal_CNV[,-1] %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(sample = rownames(.), pat = str_sub(rownames(.),1,12)) %>% 
    left_join(tcga_pscore[,c(1,4)], by = 'pat') %>% 
    filter(!is.na(Group)) %>% 
    filter(str_sub(sample,14,16) == '01A') %>% 
    filter(duplicated(.$pat) == 'FALSE') %>% 
    select(pat,sample,Group,everything()) %>% 
    pivot_longer(c(4:ncol(.)),names_to = 'Gene_Symbol', values_to = 'type') %>% 
    select(Gene_Symbol,type) %>% 
    group_by(Gene_Symbol,type) %>% 
    summarise(., count = n()) %>% 
    as_tibble() %>% 
    pivot_wider(names_from = type, values_from = count) %>% 
    mutate(`CNVs_percent  (%)` = ((as.numeric(AMP) + as.numeric(DEL)) / 
                         (as.numeric(AMP) + as.numeric(DEL) + as.numeric(None)))*100) %>% 
    select(Gene_Symbol, `CNVs_percent  (%)`)
  return(focal_CNV2)
}

focal_CNV_per <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/focal_data_by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  cnvs_per()

fabcolors = c('a' = "#D53E4F", 'b' = "#F46D43")


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F4_E.pdf',width=10,height=5)
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


save(F4A_Immune_Subtype_Box_plot,Immune_Subtype,my_comparisons,
     F4B_RNASeq_FPKM_Plot,RNASeq_FPKM,
     F4C1_Immune_Sankey_plot,Sankey_Immune,
     F4C2_IHC_Sankey_plot,Sankey_IHC,
     F4D_Lollipop_plot,wilcox_result,
     brca_maf,focal_CNV_per,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F4A_D.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F4A_D.rdata')











high_score <- as.character(subset(brca_clin, Group=="High_Score")$Tumor_Sample_Barcode)
low_score <- as.character(subset(brca_clin, Group=="Low_Score")$Tumor_Sample_Barcode)

#整理分组
high_scoreMAF <- subsetMaf(maf = allmaf, tsb = high_score)
low_scoreMAF <- subsetMaf(maf = allmaf, tsb = low_score)

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_B.pdf',width=30,height=15)
coOncoplot(m1 = high_scoreMAF, m2 = low_scoreMAF, 
          m1Name = 'High_Score', m2Name = 'Low_Score',)
dev.off()


#使用mafCompare比较差异突变基因
fvsm <- mafCompare(m1 = high_scoreMAF, m2 = low_scoreMAF, 
                   m1Name = "High_Score", m2Name = "Low_Score", minMut = 5)


#比较之后的adjPval没有达标的


#####高分低分两组CNV的差异，并绘制瀑布图#######################################################
##先整理注释文件和segments文件
if(F){
  hg_marker_file <- read.delim("/Volumes/DATA/Annotation/TCGA/snp6.na35.remap.hg38.subset.txt.gz") %>% 
    #注意：只有下载Masked Copy Number Segment是需要参数设置
    filter(freqcnv =="FALSE") %>% 
    select(c(1,2,3))
  write.table(hg_marker_file,"/Volumes/DATA/Annotation/TCGA/hg_marker_file.txt",sep = "\t",
              col.names = TRUE,row.names = F)
  
  
  load(file = "/Volumes/DATA/TCGA/TCGAbiolinks/BRCA_CNV.rda")
  tumor_seg <- eval(parse(text = load("/Volumes/DATA/TCGA/TCGAbiolinks/BRCA_CNV.rda"))) %>% 
    select('Sample','Chromosome','Start','End','Num_Probes','Segment_Mean') %>% 
    filter(substr(.$Sample,14,15) == "01") %>% 
    mutate(pat = str_sub(.$Sample, start = 1,end = 12)) %>% 
    left_join(tcga_pscore[,c(1,4)],by = 'pat')
  
  tumor_seg_high <- tumor_seg %>% 
    filter(Group == 'High_Score') %>% 
    select(-c('Group','pat'))
  
  tumor_seg_low <- tumor_seg %>% 
    filter(Group == 'Low_Score') %>% 
    select(-c('Group','pat'))
  
  write.table(tumor_seg_high,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/BRCA_CNV_high.txt',sep = '\t',quote = F,row.names = F)
  write.table(tumor_seg_low,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/BRCA_CNV_low.txt',sep = '\t',quote = F,row.names = F)
  #保存后，到在线工具网页中计算
}

##载入CNVs文件并分析，G-Score是不具备跨样本比较能力的

gistic.all <- readGistic(gisticAllLesionsFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/amp_genes.conf_90.txt", 
                         gisticDelGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/del_genes.conf_90.txt", 
                         gisticScoresFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/scores.gistic", isTCGA=TRUE)



gistic.low <- readGistic(gisticAllLesionsFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/Low_Group/all_lesions.conf_90.txt", 
                            gisticAmpGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/Low_Group/amp_genes.conf_90.txt", 
                            gisticDelGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/Low_Group/del_genes.conf_90.txt", 
                            gisticScoresFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/Low_Group/scores.gistic", isTCGA=TRUE)

gistic.high <- readGistic(gisticAllLesionsFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/High_Group/all_lesions.conf_90.txt", 
                         gisticAmpGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/High_Group/amp_genes.conf_90.txt", 
                         gisticDelGenesFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/High_Group/del_genes.conf_90.txt", 
                         gisticScoresFile="/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/High_Group/scores.gistic", isTCGA=TRUE)

gisticChromPlot(gistic = gistic.all,y_lims = c(2,-0.5))

gisticChromPlot(gistic = gistic.low,y_lims = c(2,-0.5))

gisticChromPlot(gistic = gistic.high,y_lims = c(2,-0.5))


##对CNVs进行差异分析

#整理出不同分组下的基因CNVs
##比较的内容：

#1、不同分组总数的差异
exnumber <- function(broad_CNVs) {
  broad_CNVs <- cbind(broad_CNVs[,1], 
                      as.data.frame(apply(broad_CNVs[,-1],c(1,2),function(x){
                        x = as.numeric(x)
                        #以0.3为cut-off值
                        if(x < -0.3){
                          x = 'DEL'
                        }else if(x > 0.3){
                          x = 'AMP'
                        }else{x = 'None'}
                        x = as.character(x)
                        return(x)
                      })))
  colnames(broad_CNVs) <- str_replace_all(colnames(broad_CNVs),pattern = "[.]",replacement = "-")
  rownames(broad_CNVs) <- broad_CNVs[,1]
  
  broad_CNVs <- broad_CNVs[,-1] %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(sample = rownames(.), pat = str_sub(rownames(.),1,12)) %>% 
    left_join(tcga_pscore[,c(1,4)], by = 'pat') %>% 
    filter(!is.na(Group)) %>% 
    filter(str_sub(sample,14,16) == '01A') %>% 
    filter(duplicated(.$pat) == 'FALSE') %>% 
    select(pat,sample,Group,everything()) %>% 
    pivot_longer(c(4:ncol(.)),names_to = 'Gene_Symbol', values_to = 'type') %>% 
    select(Group,Gene_Symbol,type) %>% 
    group_by(Group,Gene_Symbol,type) %>% 
    summarise(., count = n()) %>% 
    as_tibble()
}

broad_CNVs_Number <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/broad_data_by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  exnumber()

focal_CNVs_Number <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/focal_data_by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  exnumber()

all_CNVs_Number <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/all_thresholded.by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  exnumber()


#2、卡方检验找到不同分组CNVs/None的差异基因

CNV_chi.test <- function(broad_CNVs) {
  #broad_CNVs = focal_CNV_chi.test
  broad_CNVs <- cbind(broad_CNVs[,1], 
                      as.data.frame(apply(broad_CNVs[,-1],c(1,2),function(x){
                        x = as.numeric(x)
                        #以0.3为cut-off值
                        if(x < -0.3){
                          x = 'DEL'
                        }else if(x > 0.3){
                          x = 'AMP'
                        }else{x = 'None'}
                        x = as.character(x)
                        return(x)
                      })))
  
  colnames(broad_CNVs) <- str_replace_all(colnames(broad_CNVs),pattern = "[.]",replacement = "-")
  rownames(broad_CNVs) <- broad_CNVs[,1]
  
  broad_CNVs <- broad_CNVs[,-1] %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(sample = rownames(.), pat = str_sub(rownames(.),1,12)) %>% 
    left_join(tcga_pscore[,c(1,4)], by = 'pat') %>% 
    filter(!is.na(Group)) %>% 
    filter(str_sub(sample,14,16) == '01A') %>% 
    filter(duplicated(.$pat) == 'FALSE') %>% 
    select(pat,sample,Group,everything())
  
  high_score <- broad_CNVs %>% 
    filter(Group == 'High_Score') %>% 
    select(-c('pat','sample','Group')) %>% 
    apply(.,2,function(x){
      #x=high_score[,12]
      Wild = length(which(x == 'None'))
      DEL = length(which(x == 'DEL'))
      AMP = length(which(x == 'AMP'))
      CNVs = DEL + AMP
      a = rbind(CNVs,Wild)
      colnames(a) = names(x)
      
      return(a)
    })
  rownames(high_score) = c('High_CNVs','High_Wild')
  
  high_score <- high_score %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(gene = rownames(.))
  
  
  
  low_score <- broad_CNVs %>% 
    filter(Group == 'Low_Score') %>% 
    select(-c('pat','sample','Group')) %>% 
    apply(.,2,function(x){
      #x=high_score[,12]
      Wild = length(which(x == 'None'))
      DEL = length(which(x == 'DEL'))
      AMP = length(which(x == 'AMP'))
      CNVs = DEL + AMP
      a = rbind(CNVs,Wild)
      colnames(a) = names(x)
      
      return(a)
    })
  rownames(low_score) = c('Low_CNVs','Low_Wild')
  
  score <- low_score %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(gene = rownames(.)) %>% 
    left_join(high_score, by = "gene") %>% 
    select(gene,everything())
  rownames(score) <- score[,1]
  score <- score[,-1] %>% 
    as.data.frame()
  
  rt = score
  func <- function(i){
    x=matrix(c(rt[i,1],rt[i,2],rt[i,3],rt[i,4]), ncol = 2)
    chiTest=chisq.test(x)
    Low_score_Ratio=rt[i,1]/(rt[i,1]+rt[i,2])
    High_score_Ratio=rt[i,3]/(rt[i,3]+rt[i,4])
    Gene=row.names(rt[i,])
    Stat=chiTest$statistic
    Pvalue=chiTest$p.value
    outTab=cbind(Gene,Low_score_Ratio,High_score_Ratio,Stat,Pvalue)
    return(outTab)
  }
  
  
  cl <- makeCluster(6)
  registerDoParallel(cl)
  outTab <- foreach(x = 1:nrow(rt),.combine='rbind') %dopar% func(x)
  stopCluster(cl)
  
  broad_CNVs_chi <- outTab %>% 
    as.data.frame() %>% 
    mutate(adjPvalue = p.adjust(as.numeric(as.vector(.[,"Pvalue"])),method ="fdr")) %>% 
    filter(adjPvalue < 0.1)
  return(broad_CNVs_chi)
}


broad_CNV_chi.test <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/broad_data_by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  CNV_chi.test()
write_csv(broad_CNV_chi.test,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/broad_CNV_chi.test.csv')


focal_CNV_chi.test <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/focal_data_by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  CNV_chi.test()
write_csv(focal_CNV_chi.test,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/focal_CNV_chi.test.csv')


all_CNV_chi.test <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/input/GISTIC/ALL/all_thresholded.by_genes.txt") %>% 
  select(-c(2:3)) %>% 
  CNV_chi.test()
write_csv(all_CNV_chi.test,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/all_CNV_chi.test.csv')

write_csv(all_CNV_chi.test,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/7_all_CNV_chi.test.csv')




