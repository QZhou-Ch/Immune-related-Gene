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

###1 执行ESTIMATE#############

estimate <- function(dat,idtype,platform){
  #data是行名为基因名，列名为样品名
  #idtype %in% c("GeneSymbol", "EntrezID")
  #platform %in% c("affymetrix", "agilent", "illumina")
  input.f = 'estimate_input.txt'
  output.f = 'estimate_gene.gct'
  output.ds = 'estimate_score.gct'
  write.table(dat,file = input.f,sep = '\t',quote = F)
  library(estimate)
  filterCommonGenes(input.f = input.f,
                    output.f = output.f ,
                    id = idtype)
  estimateScore(input.ds = output.f,
                output.ds = output.ds,
                platform = platform)   ## 注意这个platform参数
  scores=read.table(output.ds,skip = 2,header = T)
  rownames(scores)=scores[,1]
  scores=t(scores[,3:ncol(scores)])
  return(scores)
}

estimate_GPL570 <- estimate(exp_mRNA[,-1], idtype = "EntrezID", platform = "affymetrix") %>% 
  as.data.frame() %>% 
  mutate(pat = rownames(.)) %>% 
  select('pat',everything())

write.table(estimate_GPL570, file = './table/3.2estimate_GPL570.txt',row.names = F, sep = '\t',quote = F)
#estimate_GPL570 <- read.delim("./table/3.2estimate_GPL570.txt")



#汇总数据
ALL_plot <- estimate_GPL570 %>% 
  left_join(rbind(GPL570_Train_clust,GPL570_Test_clust)[,c(1,2)], by = 'pat') %>% 
  filter(!is.na(score)) %>% 
  mutate(Group = ifelse(score < median(score),"Low_Score","High_Score")) %>% 
  select("pat","Group","score","StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity",everything())

#wilcox检验及相关分析

wilcox_result <- map(ALL_plot[,-c(1,2,3)],function(x){
  #x = ALL_plot[,-c(1,2)][,18]
  wilcox_res = wilcox.test(x~ALL_plot$Group)
  a = wilcox_res[["p.value"]]
  b = wilcox_res[["statistic"]]
  c = cbind(x,ALL_plot$Group) %>% 
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

rownames(wilcox_result) <- colnames(ALL_plot[,-c(1,2,3)])


GPL570_estimate <- cbind(ALL_plot[,c(1,2)],
                               apply(ALL_plot[,c(4:7)],2,scale)) %>% 
  pivot_longer(-c(pat,Group),names_to = "ESTIMATE", values_to = "value") %>% 
  mutate(value = as.numeric(value))


F2E_GPL570_estimate_Plot <- ggplot(GPL570_estimate, aes(x = Group, y = value,fill = Group))+
  geom_boxplot()+
  stat_compare_means(method = 'wilcox.test',size = 3.2, vjust = -0.5,hjust = -0.5)+
  facet_grid(.~ESTIMATE)+
  labs(y = "ESTIMATE", x = "")+
  theme(panel.background=element_rect(colour="black",fill=NA),
        panel.grid=element_blank(),
        legend.position="none",
        legend.background=element_rect(colour=NA,fill=NA),
        strip.text = element_text(size = 15)
  )


F2E_GPL570_estimate_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F2_D.pdf',width=30,height=15)
F2E_GPL570_estimate_Plot
dev.off()
save(F2E_GPL570_estimate_Plot,GPL570_estimate, 
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2E.rdata')

load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F2E.rdata')










