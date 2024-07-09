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

rootdir <- ""

# 1 整理免疫细胞及癌细胞的芯片数据 -------------------------------------------------------


#加载免疫细胞及癌细胞的芯片row数据
immuneLncRNA_gpl570_eset.rma <- justRMA(filenames = 
                                          dir(path = str_c(rootdir,"/DATABASE/GEO/rawdata/ICB_in_BRCA/gpl570_in_immunecell/cel")), 
                                        celfile.path=str_c(rootdir,"/DATABASE/GEO/rawdata/ICB_in_BRCA/gpl570_in_immunecell/cel"))
#整理出ccle的表达矩阵
datExpr.ccle <- exprs(immuneLncRNA_gpl570_eset.rma) %>% 
  as.data.frame() %>% 
  select(-which(str_detect(colnames(.),'GSM')))
colnames(datExpr.ccle) <- str_split(colnames(datExpr.ccle),pattern = '.cel',simplify = T)[,1]


#整理所有细胞的pdata
pdata.immunecell <- read.delim(str_c(rootdir,"/DATABASE/GEO/rawdata/ICB_in_BRCA/gpl570_in_immunecell/normal.txt"))
pdata.ccle <- cbind('ccle',colnames(datExpr.ccle),'breast')

colnames(pdata.ccle) <- colnames(pdata.immunecell)
pdata <- rbind(pdata.immunecell,pdata.ccle)
rownames(pdata) <- pdata$GSM


#整理出immune的表达矩阵
datExpr.immune <- exprs(immuneLncRNA_gpl570_eset.rma) %>% 
  as.data.frame() %>% 
  select(which(str_detect(colnames(.),'GSM')))
colnames(datExpr.immune) <- str_split(colnames(datExpr.immune),pattern = '[.]',simplify = T)[,1]
colnames(datExpr.immune) <- str_split(colnames(datExpr.immune),pattern = '[_]',simplify = T)[,1]
datExpr.immune <- datExpr.immune[,which(colnames(datExpr.immune) %in% pdata$GSM)]

#聚合两个表达矩阵
datExpr <- cbind(datExpr.ccle,datExpr.immune)%>% 
  as.data.frame() %>% 
  mutate(Probe_ID = rownames(.))

save(datExpr,datExpr.ccle,datExpr.immune,pdata,file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/immune-brca_cellline.rdata'))
#load(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/immune-brca_cellline.rdata'))



# 2 分别在癌细胞和免疫细胞中筛出GPL570平台中的lncRNA ---------------------------------------


Probe_Ensembl <- read.csv(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/probe2ensembl.csv'))
#load(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/MicroArray!!!/rdata/0gpl57_annotation.rdata'))

datExpr_lnc <-  Probe_Ensembl %>% 
  filter(gene_biotype == 'lncRNA') %>% 
  left_join(datExpr, by = 'Probe_ID') %>% 
  select(-'gene_biotype') %>% 
  as.data.frame() 

#immune lncRNA expression
datExpr.immune.lnc <- datExpr_lnc[,c(1,2,which(colnames(datExpr_lnc) %in% colnames(datExpr.immune)))]

#一对多的情况下，取最大值。
datExpr.immune.lnc <- aggregate(x = datExpr.immune.lnc[,3:ncol(datExpr.immune.lnc)],
                                by = list(datExpr.immune.lnc$Ensembl),
                                FUN = max) %>% 
  column_to_rownames(., var = "Group.1")


#ccle lncRNA expression
datExpr.ccle.lnc <- datExpr_lnc[,c(1,2,which(colnames(datExpr_lnc) %in% colnames(datExpr.ccle)))]

#一对多的情况下，取最大值。
datExpr.ccle.lnc <- aggregate(x = datExpr.ccle.lnc[,3:ncol(datExpr.ccle.lnc)],
                                by = list(datExpr.ccle.lnc$Ensembl),
                                FUN = max) %>% 
  column_to_rownames(., var = "Group.1")

# 3 构造高表达和低表达的函数 ----------------------------------------------------------


#每种细胞中的前/后5%表达量，然后union或者intersect
highexpression <- function(expr,cutoff){
  #expr <- immunecell_miRNA
  #cutoff = 0.05
  #首先保留每种细胞前5%最高表达量
  expr <- as.data.frame(expr)
  name <- map(expr,function(x){
    expr <- expr[order(x,decreasing = T),]
    name <- rownames(expr)[1:round(nrow(expr)*cutoff)]
    return(name)
  }) %>% 
    #purrr::reduce(intersect)
    purrr::reduce(union)
  return(name)
}


lowexpression <- function(expr,cutoff){
  #expr <- immunecell_miRNA
  #cutoff = 0.05
  #首先保留每种细胞前5%最低表达量
  expr <- as.data.frame(expr)
  name <- map(expr,function(x){
    expr <- expr[order(x,decreasing = F),]
    name <- rownames(expr)[1:round(nrow(expr)*cutoff)]
    return(name)
  }) %>% 
    purrr::reduce(intersect)
    #purrr::reduce(union)
  return(name)
}


# 4 执行 --------------------------------------------------------------------


immune_lncrna <- intersect(highexpression(datExpr.immune.lnc,0.2),
                           lowexpression(datExpr.ccle.lnc,0.2))

write.csv(immune_lncrna,file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/4_immune_lncRNA.csv'))

save(immuneLncRNA_gpl570_eset.rma,pdata,immune_lncrna,Probe_Ensembl,
     file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/1.2immuneLncRNA_NEW_ARTICLE_gpl.rdata'))
  

#load(file = '/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/1.2immuneLncRNA_NEW_ARTICLE_gpl.rdata')



