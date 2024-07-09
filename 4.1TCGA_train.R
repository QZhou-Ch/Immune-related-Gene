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

load('/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/3.1GPL570_train.rdata')


#####1 依据counts整理出所需的表达矩阵和生存信息###################

if(T){
  load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/2.1TCGA_counts.rdata')
  tcga_pdata <- tcga_BR[,c(1:67)]
  
  tcga_rnaseq <- tcga_BR[,c(2,68:ncol(tcga_BR))] %>% 
    as.data.frame()
  rownames(tcga_rnaseq) <- tcga_rnaseq$sample
  
  tcga_rnaseq <- cbind(tcga_rnaseq$sample,t(limma::voom(t(as.matrix(tcga_rnaseq[,-1])))$E))
  colnames(tcga_rnaseq)[1] <- "sample"
  
  tcga_OS <- tcga_rnaseq[,-1] %>% 
    as.data.frame() %>% 
    filter(str_sub(rownames(.),16,16) == 'A') %>% 
    filter(!(rownames(.) %in% c('TCGA-A7-A0DB-01A-11R-A00Z-07','TCGA-A7-A13D-01A-13R-A12P-07',
                                'TCGA-A7-A13E-01A-11R-A12P-07','TCGA-A7-A26E-01A-11R-A169-07',
                                'TCGA-A7-A26J-01A-11R-A169-07'))) %>% 
    mutate(sample = rownames(.)) %>% 
    left_join(tcga_pdata[,c(2:4)],by = 'sample') %>% 
    filter(!(is.na(OS.Time) | is.na(OS))) %>% 
    filter(!(OS.Time == 0 | OS.Time == '#N/A')) %>% 
    #filter(!(OS == 0 & OS.Time < 365)) %>% 
    select(sample,OS,OS.Time,which(colnames(.) %in% Gene_coef$ENSEMBL))
  rownames(tcga_OS) <- tcga_OS[,1]
  
  tcga_OS <- apply(tcga_OS[,-1], c(1,2), as.numeric) %>% 
    as.data.frame()
}

#####2 计算每个sample的分数#######################################

res_OScox <- coxph(Surv(OS.Time, OS) ~ ., data =  tcga_OS)
summary(res_OScox)

res <- summary(res_OScox)[["coefficients"]] %>% 
  as.data.frame() 


OS_genename <- cbind(str_remove_all(rownames(res),pattern = '[`]'),res$coef) %>% 
  as.data.frame() %>% 
  rename(genename = V1, coef = V2)
write.csv(OS_genename, file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/6_20Gene_and_Coef_TCGA.csv')

#计算每个sample的分数
tcga_score <- tcga_OS[,-c(1,2)] %>% 
  #select(Gene_coef[,1]) %>% 
  select(OS_genename[,1]) %>% 
  apply(.,1,function(x){
  #x = tcga_OS[,-c(1,2)][1,]
  #计算两个向量的内积
  #tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(Gene_coef[,2]))
  tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(OS_genename[,2]))
  return(tt)
})


tcga_clust <- cbind(names(tcga_score),as.numeric(tcga_score)) %>% 
  as_tibble() %>% 
  rename(sample = V1,
         score = V2) %>% 
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
         time = as.numeric(tcga_OS$OS.Time),
         status = as.numeric(tcga_OS$OS)) %>% 
  mutate(score = as.numeric(score),
         time = as.numeric(time),
         status = as.numeric(status))



#####绘制K-M曲线及AUC曲线#############################################

#绘制K-M曲线
TCGA_fit_OS <- survfit(Surv(time, status) ~Group, data = tcga_clust)


logrank_OS <- survdiff(Surv(time, status) ~ Group, data = tcga_clust)
logrank_OS
pval_OS = 1 - pchisq(logrank_OS$chisq, length(logrank_OS$n) - 1)
pval_OS
HR_OS = (logrank_OS$obs[2]/logrank_OS$exp[2])/(logrank_OS$obs[1]/logrank_OS$exp[1])
HR_OS
low95_OS = exp(log(HR_OS) - qnorm(0.975)*sqrt(1/logrank_OS$exp[2]+1/logrank_OS$exp[1]))
low95_OS
up95_OS = exp(log(HR_OS) + qnorm(0.975)*sqrt(1/logrank_OS$exp[2]+1/logrank_OS$exp[1]))
up95_OS


F3A_TCGA_KM_OS_Plot <- ggsurvplot(TCGA_fit_OS, data = tcga_clust, 
           risk.table = T, conf.int = T,pval = T,
           xlab = "Follow up time(DAY)",# 指定x轴标签
           ylab = "OS Probability",
           #legend = c(0.8,0.75),# 指定图例位置
           break.x.by = 1000,  # 设置x轴刻度间距
           xlim = c(0,4000))
F3A_TCGA_KM_OS_Plot


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F3_A.pdf',width=15,height=15)
F3A_TCGA_KM_OS_Plot
dev.off()


multi_ROC <- function(time_vector,time_name,risk_score_table){
  single_ROC <- function(single_time,single_time_name){
    for_ROC <- survivalROC(Stime = risk_score_table$time,
                           status = risk_score_table$status,
                           marker = risk_score_table$score,
                           predict.time = single_time, method = 'KM')
    data.frame('Time_point'=rep(single_time_name, length(for_ROC$TP)),
               'True_positive'=for_ROC$TP,
               'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values,
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- pmap(list(time_vector, time_name), single_ROC) %>% 
    reduce(rbind) %>% 
    mutate(Time_point = factor(.$Time_point,levels = time_name))
}


time_vector <- 365*c(3,5,10)
time_name <- c('3y','5y','10y')

TCGA_OS_clust_ROC <- multi_ROC(time_vector = time_vector,time_name = time_name,risk_score_table = tcga_clust)



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
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))
F3B_TCGA_OS_ROC_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F3_B.pdf',width=15,height=15)
F3B_TCGA_OS_ROC_Plot
dev.off()


#####根据PFS绘制KM曲线和AUC曲线##########################################

tcga_PFS <- tcga_rnaseq[,-1] %>% 
  t() %>% 
  as.data.frame() %>% 
  select(which(str_sub(colnames(.),16,16) == 'A'))  %>% 
  select(-c('TCGA-A7-A0DB-01A-11R-A00Z-07','TCGA-A7-A13D-01A-13R-A12P-07',
            'TCGA-A7-A13E-01A-11R-A12P-07','TCGA-A7-A26E-01A-11R-A169-07',
            'TCGA-A7-A26J-01A-11R-A169-07')) %>%
  t() %>% 
  as.data.frame() %>% 
  mutate(sample = rownames(.)) %>% 
  left_join(tcga_pdata[,c(2,5,6)],by = 'sample') %>% 
  filter(!(is.na(PFI.Time) | is.na(PFI))) %>% 
  filter(!(PFI.Time == 0 | PFI.Time == '#N/A')) %>% 
  #filter(!(PFI == 0 & PFI.Time < 365)) %>% 
  select(sample,PFI,PFI.Time,which(colnames(.) %in% Gene_coef$ENSEMBL))
rownames(tcga_PFS) <- tcga_PFS[,1]

tcga_PFS <- apply(tcga_PFS[,-1], c(1,2), as.numeric) %>% 
  as.data.frame()


PFI_score <- tcga_score


PFI_clust <- cbind(names(PFI_score),as.numeric(PFI_score)) %>% 
  as_tibble() %>% 
  rename(sample = V1,
         score = V2) %>% 
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
         time = as.numeric(tcga_PFS$PFI.Time),
         status = as.numeric(tcga_PFS$PFI)) %>% 
  mutate(score = as.numeric(score),
         time = as.numeric(time),
         status = as.numeric(status))


TCGA_fit_PFI <- survfit(Surv(time, status) ~Group, data = PFI_clust)

logrank_PFI <- survdiff(Surv(time, status) ~ Group, data = PFI_clust)

pval_PFI = 1 - pchisq(logrank_PFI$chisq, length(logrank_PFI$n) - 1)
pval_PFI
HR_PFI = (logrank_PFI$obs[2]/logrank_PFI$exp[2])/(logrank_PFI$obs[1]/logrank_PFI$exp[1])
HR_PFI
low95_PFI = exp(log(HR_PFI) - qnorm(0.975)*sqrt(1/logrank_PFI$exp[2]+1/logrank_PFI$exp[1]))
low95_PFI
up95_PFI = exp(log(HR_PFI) + qnorm(0.975)*sqrt(1/logrank_PFI$exp[2]+1/logrank_PFI$exp[1]))
up95_PFI


F3C_TCGA_KM_PFI_Plot <- ggsurvplot(TCGA_fit_PFI, data = PFI_clust, 
           risk.table = T, conf.int = T,pval = T,
           xlab = "Follow up time(DAY)",# 指定x轴标签
           ylab = "PFS Probability",
           #legend = c(0.8,0.75),# 指定图例位置
           break.x.by = 1000,  # 设置x轴刻度间距
           xlim = c(0,4000))
F3C_TCGA_KM_PFI_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F3_C.pdf',width=15,height=15)
F3C_TCGA_KM_PFI_Plot
dev.off()



TCGA_PFI_clust_ROC <- multi_ROC(time_vector = time_vector,time_name = time_name,risk_score_table = PFI_clust)



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
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))
F3D_TCGA_PFI_ROC_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F3_D.pdf',width=15,height=15)
F3D_TCGA_PFI_ROC_Plot
dev.off()



#####计算所有患者的分数################################################


tcga_pscore <- tcga_rnaseq[,which(colnames(tcga_rnaseq) %in% OS_genename$genename)] %>% 
  t() %>% 
  as.data.frame() %>% 
  select(which(str_sub(colnames(.),16,16) == 'A'))  %>% 
  select(-c('TCGA-A7-A0DB-01A-11R-A00Z-07','TCGA-A7-A13D-01A-13R-A12P-07',
            'TCGA-A7-A13E-01A-11R-A12P-07','TCGA-A7-A26E-01A-11R-A169-07',
            'TCGA-A7-A26J-01A-11R-A169-07')) %>%
  as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  select(OS_genename[,1]) %>% 
  apply(.,1,function(x){
    tt = as.numeric(as_vector(x)) %*% as.numeric(as_vector(OS_genename[,2]))
    return(tt)
  })

tcga_pscore <- cbind(names(tcga_pscore),as.numeric(tcga_pscore)) %>% 
  as_tibble() %>% 
  rename(sample = V1,
         score = V2) %>% 
  left_join(tcga_pdata, by = 'sample') %>% 
  select(colnames(tcga_pdata),everything()) %>% 
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score")) %>% 
  select(pat,sample,score,Group,everything())


#####整理出分组患者及分数########################################################

load(file = '/Volumes/DATA/rdata/ICB_in_BRCA/RNASeq/2.1TCGA_FPKM.rdata')
tcga_pdata <- tcga_BR[,c(1:67)]

tcga_rnaseq_FPKM<- tcga_BR[,c(2,68:ncol(tcga_BR))] %>% 
  as.data.frame()
rownames(tcga_rnaseq_FPKM) <- tcga_rnaseq_FPKM$sample



#####评分，年龄，stage，与OS，PFI做多因素生存分析，并绘制森林图################################

TCGA_multiCOX_OS <- tcga_pscore[,c(1,3,5:10)] %>% 
  filter(!(stage == "NA")) %>% 
  filter(!(is.na(OS.Time) | is.na(OS))) %>% 
  filter(!(OS.Time == 0 | OS.Time == '#N/A')) %>% 
  #filter(!(OS == 0 & OS.Time < 365)) %>% 
  mutate(#early_stage = ifelse(stage == '1',"I","II/III"),
    #meta_stage = ifelse(stage == '4',"IV","I/II/III"),
    OS.Time = as.numeric(OS.Time),
    score = as.numeric(score),
    age = as.numeric(age),
    stage = factor(stage,labels = c('I','II','III','IV')),
    OS = as.numeric(OS)) %>% 
  select(OS.Time,OS,age,stage,score)

TCGA_OS_model <- coxph( Surv(OS.Time, OS) ~  age + stage + score , data =  TCGA_multiCOX_OS)
summary(TCGA_OS_model)
TCGA_OS_model


F3E_TCGA_OS_forest_Plot <- ggforest(TCGA_OS_model,  #coxph得到的Cox回归结果
         data = TCGA_multiCOX_OS,  #数据集
         #标题
         main = 'Hazard ratio of BRCA (OS)',  
         #前三列距离
         cpositions = c(0.05, 0.15, 0.35),  
         #字体大小
         fontsize = 1, 
         #相对变量的数值标签，也可改为1
         refLabel = 'reference', 
         #保留HR值以及95%CI的小数位数
         noDigits = 3)+
  theme_void()+
  theme(strip.background = element_rect(fill="grey90"),
        strip.text = element_text(size=13,face="plain",color="black"),
        axis.title=element_text(size=13,face="plain",color="black"),
        axis.text = element_text(size=11,face="plain",color="black"),
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        axis.ticks=element_line(colour="black"))
F3E_TCGA_OS_forest_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_I.pdf',width=30,height=15)
F3E_TCGA_OS_forest_Plot
dev.off()




TCGA_multiCOX_PFI <- tcga_pscore[,c(1,3,7:10)] %>% 
  filter(!(stage == "NA")) %>% 
  filter(!(is.na(PFI.Time) | is.na(PFI))) %>% 
  filter(!(PFI.Time == 0 | PFI.Time == '#N/A')) %>% 
  #filter(!(OS == 0 & OS.Time < 365)) %>% 
  mutate(#early_stage = ifelse(stage == '1',"I","II/III"),
    #meta_stage = ifelse(stage == '4',"IV","I/II/III"),
    PFI.Time = as.numeric(PFI.Time),
    score = as.numeric(score),
    age = as.numeric(age),
    stage = factor(stage,labels = c('I','II','III','IV')),
    PFI = as.numeric(PFI)) %>% 
  select(PFI.Time,PFI,age,stage,score)

TCGA_PFI_model <- coxph( Surv(PFI.Time, PFI) ~  age + stage + score , data =  TCGA_multiCOX_PFI)
summary(TCGA_PFI_model)
TCGA_PFI_model

F3F_TCGA_PFI_forest_Plot <- ggforest(TCGA_PFI_model,  #coxph得到的Cox回归结果
         data = TCGA_multiCOX_PFI,  #数据集
         #标题
         main = 'Hazard ratio of BRCA (PFS)',  
         #前三列距离
         cpositions = c(0.05, 0.15, 0.35),  
         #字体大小
         fontsize = 1, 
         #相对变量的数值标签，也可改为1
         refLabel = 'reference', 
         #保留HR值以及95%CI的小数位数
         noDigits = 3)+
  theme_void()+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA))
F3F_TCGA_PFI_forest_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F3_F.pdf',width=30,height=15)
F3F_TCGA_PFI_forest_Plot
dev.off()



save(tcga_OS,
     tcga_pscore,tcga_rnaseq_FPKM,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/4.1TCGA_Score.rdata')
#load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/4.1TCGA_Score.rdata')

save(F3A_TCGA_KM_OS_Plot,TCGA_fit_OS,tcga_clust, 
     F3B_TCGA_OS_ROC_Plot,TCGA_OS_clust_ROC,time_name,
     F3C_TCGA_KM_PFI_Plot,TCGA_fit_PFI,PFI_clust, 
     F3D_TCGA_PFI_ROC_Plot,TCGA_PFI_clust_ROC,
     F3E_TCGA_OS_forest_Plot,TCGA_OS_model, TCGA_multiCOX_OS,
     F3F_TCGA_PFI_forest_Plot,TCGA_PFI_model, TCGA_multiCOX_PFI,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F3A_F.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F3A_F.rdata')


