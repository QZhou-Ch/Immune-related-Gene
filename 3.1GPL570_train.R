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

load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/3.0renew.rdata')#整理出所需的表达矩阵和生存信息

if(T){
  #整理出lnc和mRNA
  #####0 整理出GPL570乳腺癌数据集#######################
  GPL570BR_pdata %<>% 
    filter(!(is.na(OS.Time) | is.na(OS))) %<>% 
    filter(!(OS.Time == "NA" | OS == "NA"))
write.csv(GPL570BR_pdata,file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/9_GPL570_pDATA.csv')
  
  #####1 在训练集中整理出生存相关的lncRNA#######################
  GPL570_lncRNA_sur <- exp_lncRNA %>% 
    t()
  colnames(GPL570_lncRNA_sur) <- GPL570_lncRNA_sur[1,]
  
  GPL570_lncRNA_sur <- GPL570_lncRNA_sur[-1,] %>% 
    as.data.frame() %>% 
    mutate(pat = rownames(.))
  
  GPL570_lncRNA_sur <- GPL570BR_pdata[,c(1,10,11)] %>% 
    left_join(GPL570_lncRNA_sur, by = 'pat')
  
  rownames(GPL570_lncRNA_sur) <- GPL570_lncRNA_sur$pat
  
  
  #使用单因素COX整理生存相关lncRNA
  data <- apply(GPL570_lncRNA_sur[,-1], c(1,2), as.numeric) %>% 
    as.data.frame()
  
  table(GPL570BR_pdata$batch)
  #整理出训练组
  data <- data[which(rownames(data) %in%
                       rownames(GPL570BR_pdata)
                     [-which(GPL570BR_pdata$batch %in% c("GSE20685"))]),]
  
  
  #构建函数:需传入两个变量，x是单因素，
  #expr是表达矩阵，其中第二列是生存时间，第一列是状态
  
  coxf<-function(x,expr){
    fmla1 <- as.formula(Surv(expr[,2],expr[,1])~expr[,x])
    mycox <- coxph(fmla1,data = expr)
  }
  
  
  #单因素cox回归
  #使用for循环做单因素DFS-Cox生存分析
  coxR <- data.frame()
  for(a in colnames(data[,3:ncol(data)])){
    mycox=coxf(a,data)
    coxResult = summary(mycox)
    coxR=rbind(coxR,cbind(at.name=a,HR=coxResult$coefficients[,"exp(coef)"],
                          P=coxResult$coefficients[,"Pr(>|z|)"]))
  }
  
  
  GPL570_lncRNA <- exp_lncRNA[which(exp_lncRNA$Probe_ID %in% 
                                      (coxR[which(coxR$P < 0.05),1])),] %>% 
    rename(genename = Probe_ID) %>% 
    select(genename,everything())
  rownames(GPL570_lncRNA) <- GPL570_lncRNA$genename
  
  
  #####2 在训练集中整理出生存相关的mRNA#######################
  
  GPL570_mRNA_sur <- exp_mRNA[which(exp_mRNA$id %in% ibc_therapy_mRNA$ENTREZID),] %>% 
    as.data.frame() %>% 
    t()
  colnames(GPL570_mRNA_sur) <- GPL570_mRNA_sur[1,]
  
  GPL570_mRNA_sur <- GPL570_mRNA_sur[-1,] %>% 
    as.data.frame() %>% 
    mutate(pat = rownames(.))
  
  GPL570_mRNA_sur <- GPL570BR_pdata[,c(1,10,11)] %>% 
    left_join(GPL570_mRNA_sur, by = 'pat')
  
  rownames(GPL570_mRNA_sur) <- GPL570_mRNA_sur$pat
  
  
  #使用单因素COX整理生存相关mRNA
  data <- apply(GPL570_mRNA_sur[,-1], c(1,2), as.numeric) %>% 
    as.data.frame()
  #整理出训练组
  data <- data[which(rownames(data) %in%
                       rownames(GPL570BR_pdata)
                     [-which(GPL570BR_pdata$batch %in% c("GSE20685"))]),]
  
  
  #单因素cox回归
  #使用for循环做单因素OS-Cox生存分析
  coxR <- data.frame()
  for(a in colnames(data[,3:ncol(data)])){
    mycox=coxf(a,data)
    coxResult = summary(mycox)
    coxR=rbind(coxR,cbind(at.name=a,HR=coxResult$coefficients[,"exp(coef)"],
                          P=coxResult$coefficients[,"Pr(>|z|)"]))
  }
  
  
  
  GPL570_mRNA <- exp_mRNA[which(exp_mRNA$id %in% 
                                  (coxR[which(coxR$P < 0.05),1])),] %>% 
    left_join(clusterProfiler::bitr(.$id,fromType = "ENTREZID", 
                                    toType = "SYMBOL", 
                                    OrgDb = "org.Hs.eg.db"),
              by = c('id' = 'ENTREZID')) %>% 
    rename(genename = SYMBOL) %>% 
    select(-id) %>% 
    select(genename,everything()) #%>% filter(!(genename %in% c('TEK','EGFR')))
  
  rownames(GPL570_mRNA) <- GPL570_mRNA$genename
  
  
  
  #####3 整合mRNA和lncRNA准备进行下一步的生存分析#######################
  
  #####3.1 整理数据集
  
  GPL570_mCOX <- rbind(GPL570_mRNA,GPL570_lncRNA)[,-1] %>%
    as.matrix() %>% 
    apply(.,c(1,2),as.numeric) %>% 
    t() %>% 
    as.data.frame() %>% 
    mutate(pat = rownames(.)) %>% 
    left_join(GPL570BR_pdata[,c(1,10,11)],by = 'pat') %>% 
    select(pat,OS,OS.Time,everything())
  rownames(GPL570_mCOX) <- GPL570_mCOX[,1]
  
  GPL570_mCOX <- apply(GPL570_mCOX[,-1], c(1,2), as.numeric) %>% 
    as.data.frame() %>% 
    filter(!(is.na(OS.Time) | is.na(OS))) %>% 
    #filter(!(OS == 0 & OS.Time < 365)) %>% 
    filter(!(OS.Time == 0))
  
  #####3.2 筛选lncRNA
  #####对lncRNA及mRNA进行相关分析检验，排除存在相关关系的lncRNA
  #首先就是接受正态分布shapiro.test检验
  shapiro_test <- apply(GPL570_mCOX[,-c(1,2)],2,function(x){
    pvalue <- shapiro.test(x)[["p.value"]]
    return(pvalue)
  }) %>% 
    .[which(.>0.05)]
  
  #进行spearman的相关分析
  cor_lnc2m <- GPL570_lncRNA$genename %>% 
    map(function(x){
      #x = GPL570_lncRNA$genename[1]
      lncname <- x 
      cor_spreaman <- GPL570_mCOX[,which(colnames(GPL570_mCOX) %in% GPL570_mRNA$genename)] %>% 
        map(function(x){
          a <- cor.test(x,GPL570_mCOX[,lncname], method = "spearman")
          matrix <- c(a[["p.value"]],a[["estimate"]][["rho"]])
          return(matrix)
        }) %>% 
        purrr::reduce(rbind) %>% 
        as.data.frame()
      colnames(cor_spreaman) <- c('pvalue','rho')
      rownames(cor_spreaman) <- colnames(GPL570_mCOX[,which(colnames(GPL570_mCOX) %in% GPL570_mRNA$genename)])
      #选择rho>0.5
      cor_spreaman <- cor_spreaman[which(abs(cor_spreaman$rho)>0.5 & cor_spreaman$pvalue < 0.05),]
      return(cor_spreaman)
    })
  names(cor_lnc2m) <- GPL570_lncRNA$genename
  #根据cor_lnc2m的结果，去除rho>0.5，pvalue <0.05的lncRNA。
  nocor_lnc <- cor_lnc2m %>% 
    map2(names(.),function(x,y){
      if(length(rownames(x)) == 0){
        return(y)
      }
    }) %>% 
    purrr::reduce(c)
  
  GPL570_mCOX <- GPL570_mCOX[,c(1,2,which(colnames(GPL570_mCOX) %in% GPL570_mRNA$genename),
                                which(colnames(GPL570_mCOX) %in% nocor_lnc))] %>% 
    mutate(pat = rownames(.)) %>% 
    left_join(GPL570BR_pdata[,c(1,12)], by = 'pat') %>% 
    select('pat','batch',everything())
  rownames(GPL570_mCOX) <- GPL570_mCOX$pat
  GPL570_mCOX <- GPL570_mCOX[,-1]
  
  gpl570lnc_Annotation[which(gpl570lnc_Annotation$Probe_ID %in% colnames(GPL570_mCOX)),]
  
  
  #手动删除不需要的lncRNA
  GPL570_mCOX <- GPL570_mCOX[,-which(colnames(GPL570_mCOX) %in% c("237075_at"))]
  
}


#####4 在训练集中训练，评估参数，得到Immune_Score#######################
table(GPL570_mCOX$batch)

#####4.1 建立训练集，在训练集上针对筛选RNA进行mutiCOX并计算Immune_Score
GPL570_mCOX.Train <- GPL570_mCOX[-which(GPL570BR_pdata$batch %in% c("GSE20685")),] %>% 
  select(-'batch')

res.cox <- coxph(Surv(OS.Time, OS) ~ ., data =  GPL570_mCOX.Train)
summary(res.cox)

res <- summary(res.cox)[["coefficients"]] %>% 
  as.data.frame() #%>% filter(`Pr(>|z|)` < 0.5)

mCOX_genename <- cbind(str_remove_all(rownames(res),pattern = '[`]'),res$coef) %>% 
  as.data.frame() %>% 
  rename(genename = V1, coef = V2)

#计算每个sample的分数
GPL570_score <- GPL570_mCOX.Train[,-c(1,2)] %>% 
  select(mCOX_genename[,1]) %>% 
  apply(.,1,function(x){
    #x = tcga_lasoCOX[,-c(1,2)][1,]
    #计算两个向量的内积
    tt = as_vector(x) %*% as.numeric(as_vector(mCOX_genename[,2]))
    return(tt)
  })


#####4.2 在训练集上，检验Immune_Score的预测效果


#在训练集上根据Immune_Score的中位数分组
GPL570_Train_clust <- cbind(names(GPL570_score),as.numeric(GPL570_score)) %>%
  as_tibble() %>%
  rename(pat = V1, score = V2) %>%
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
         time = as.numeric(GPL570_mCOX.Train$OS.Time),
         status = as.numeric(GPL570_mCOX.Train$OS)) %>% 
  mutate(score = as.numeric(score),
         time = as.numeric(time),
         status = as.numeric(status))

#绘制K-M曲线
GPL570_KMfit_Train <- survfit(Surv(time, status) ~Group, data = GPL570_Train_clust)
summary(GPL570_KMfit_Train)

logrank_Train <- survdiff(Surv(time, status) ~ Group, data = GPL570_Train_clust)
summary(logrank_Train)

#结果：在训练集上评估参数后的K-M生存检验
pval_Train = 1 - pchisq(logrank_Train$chisq, length(logrank_Train$n) - 1)
pval_Train
HR_Train = (logrank_Train$obs[2]/logrank_Train$exp[2])/(logrank_Train$obs[1]/logrank_Train$exp[1])
HR_Train
low95_Train = exp(log(HR_Train) - qnorm(0.975)*sqrt(1/logrank_Train$exp[2]+1/logrank_Train$exp[1]))
low95_Train
up95_Train = exp(log(HR_Train) + qnorm(0.975)*sqrt(1/logrank_Train$exp[2]+1/logrank_Train$exp[1]))
up95_Train

#绘制F2A
F2A_GPL570_KM_Train_Plot <- ggsurvplot(GPL570_KMfit_Train, data = GPL570_Train_clust, 
                                   risk.table = T, conf.int = T,pval = T,
                                   #surv.median.line = "hv",#添加中位生存线
                                   xlab = "Follow up time(DAY)",# 指定x轴标签
                                   ylab = "OS Probability",
                                   #legend = c(0.8,0.75),# 指定图例位置
                                   xlim = c(0,4000),
                                   break.x.by = 1000)
F2A_GPL570_KM_Train_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_A.pdf',width=15,height=15)
F2A_GPL570_KM_Train_Plot
dev.off()


#绘制AUC
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


GPL570_Train_clust_ROC <- multi_ROC(time_vector = time_vector,time_name = time_name,risk_score_table = GPL570_Train_clust)


###F2B
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
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))

F2B_GPL570_Train_ROC_Plot



pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_B.pdf',width=15,height=15)
F2B_GPL570_Train_ROC_Plot
dev.off()


#####5 测试集上检验Immune_Score的Rubost#######################

GPL570_mCOX.Test <- GPL570_mCOX[which(GPL570BR_pdata$batch %in% c("GSE20685")),] %>% 
  select(-'batch')

#使用训练集系数计算每个sample的分数
GPL570_score <- GPL570_mCOX.Test[,-c(1,2)] %>% 
  select(mCOX_genename[,1]) %>% 
  apply(.,1,function(x){
    #x = tcga_lasoCOX[,-c(1,2)][1,]
    #计算两个向量的内积
    tt = as_vector(x) %*% as.numeric(as_vector(mCOX_genename[,2]))
    return(tt)
  })


GPL570_Test_clust <- cbind(names(GPL570_score),as.numeric(GPL570_score)) %>%
  as_tibble() %>%
  rename(pat = V1, score = V2) %>%
  mutate(Group = ifelse(as.numeric(score) <= median(as.numeric(score)),"Low_Score","High_Score"),
         time = as.numeric(GPL570_mCOX.Test$OS.Time),
         status = as.numeric(GPL570_mCOX.Test$OS)) %>% 
  mutate(score = as.numeric(score),
         time = as.numeric(time),
         status = as.numeric(status))

#绘制K-M曲线
GPL570_KMfit_Test <- survfit(Surv(time, status) ~Group, data = GPL570_Test_clust)
GPL570_KMfit_Test
summary(GPL570_KMfit_Test)

logrank_Test <- survdiff(Surv(time, status) ~ Group, data = GPL570_Test_clust)
logrank_Test


pval_Test = 1 - pchisq(logrank_Test$chisq, length(logrank_Test$n) - 1)
pval_Test
HR_Test = (logrank_Test$obs[2]/logrank_Test$exp[2])/(logrank_Test$obs[1]/logrank_Test$exp[1])
HR_Test
low95_Test = exp(log(HR_Test) - qnorm(0.975)*sqrt(1/logrank_Test$exp[2]+1/logrank_Test$exp[1]))
low95_Test
up95_Test = exp(log(HR_Test) + qnorm(0.975)*sqrt(1/logrank_Test$exp[2]+1/logrank_Test$exp[1]))
up95_Test



##F2C
F2C_GPL570_KM_Test_Plot <- ggsurvplot(GPL570_KMfit_Test, data = GPL570_Test_clust, 
                                  risk.table = T, conf.int = T,pval = T,
                                  #surv.median.line = "hv",#添加中位生存线
                                  xlab = "Follow up time(DAY)",# 指定x轴标签
                                  ylab = "OS Probability",
                                  #legend = c(0.8,0.75),# 指定图例位置
                                  xlim = c(0,4000),
                                  break.x.by = 1000)
F2C_GPL570_KM_Test_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_C.pdf',width=15,height=15)
F2C_GPL570_KM_Test_Plot
dev.off()

#绘制AUC

GPL570_Test_clust_ROC <- multi_ROC(time_vector = time_vector,time_name = time_name,risk_score_table = GPL570_Test_clust)


#F2D
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
        panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.text = element_text(size = 16),
        legend.position=c(0.95,0.05),
        legend.justification = c(1,0),
        legend.background=element_blank(),
        legend.key = element_blank(),
        axis.ticks=element_line(colour="black"))

F2D_GPL570_Test_ROC_Plot


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_D.pdf',width=15,height=15)
F2D_GPL570_Test_ROC_Plot
dev.off()




#####6 保存结果#######################

mCOX_genename$ENSEMBL <- mCOX_genename$genename
#CCL5有两个ENSEMBL id，去掉一个
Gene_coef <- cbind(c(clusterProfiler::bitr(mCOX_genename$genename[1:15],fromType = "SYMBOL", 
                                           toType = "ENSEMBL", 
                                           OrgDb = "org.Hs.eg.db")[,2][-c(2:17,25)],
                     gpl570lnc_Annotation[which(gpl570lnc_Annotation$Probe_ID %in% colnames(GPL570_mCOX)),][,2]),
                   (mCOX_genename$coef),
                   c(mCOX_genename$genename[1:15],gpl570lnc_Annotation[which(gpl570lnc_Annotation$Probe_ID %in% colnames(GPL570_mCOX)),][,3])) %>% 
  as.data.frame() %>% 
  rename(ENSEMBL = V1, coef = V2, genename= V3) %>% 
  select('genename','ENSEMBL','coef')

write.csv(Gene_coef, file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/5_20Gene_and_Coef_GEO.csv')

save(Gene_coef,
     GPL570_Train_clust,GPL570_Test_clust,
     exp_mRNA,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/3.1GPL570_train.rdata')

#load('/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/3.1GPL570_train.rdata')


save(F2A_GPL570_KM_Train_Plot,GPL570_KMfit_Train,GPL570_Train_clust, 
     F2B_GPL570_Train_ROC_Plot,GPL570_Train_clust_ROC,time_name,
     F2C_GPL570_KM_Test_Plot,GPL570_KMfit_Test,GPL570_Test_clust, 
     F2D_GPL570_Test_ROC_Plot,GPL570_Test_clust_ROC,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2A_D.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2A_D.rdata')








