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


rootdir <- "/home/zhou"

load(file = str_c(rootdir,"/DATABASE/RProject/GSE_DATA/1.breast_with_survival_data.rdata"))



# 2 整合GEO表达数据 -------------------------------------------------------------


##2.1 RMA normalization

gpl570_BR_eset.rma <- justRMA(filenames = 
                                dir(path = str_c(rootdir,"/DATABASE/GEO/rawdata/ICB_in_BRCA/BR_Survival.gpl570")), 
                              celfile.path= 
                                str_c(rootdir,"/DATABASE/GEO/rawdata/ICB_in_BRCA/BR_Survival.gpl570"))

datExpr <- exprs(gpl570_BR_eset.rma)

#整理GSM的名字
name <- str_split(colnames(datExpr),pattern = "[.]",simplify = T)[,1]
name <- str_split(name,pattern = "_",simplify = T)[,1]
colnames(datExpr) <- name

#2.2结合pdata进行批间差矫正
meta_exp <- datExpr[,which(colnames(datExpr) %in% GPL570BR_pdata$pat)]
ComBat.meta_exp <- sva::ComBat(dat = meta_exp, batch = GPL570BR_pdata$batch)
#rbe.meta_exp <- limma::removeBatchEffect(meta_exp, batch = pdata$batch)

#gpl570_BR.pdata <- pdata[-which(pdata$m == 1 | pdata$OS_status == 'NA' | pdata$OS == 'NA'),c(1,10,11)]
GPL570_BR.expr <- ComBat.meta_exp

#####3 整理乳腺癌数据集中entrez gene_id的表达###########################################################
#合并///





#####4 整理乳腺癌数据集中lncRNA的表达##################################################################

#load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/0gpl57_annotation.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/1.2immuneLncRNA_microarray_gpl.rdata')


exp_lncRNA <- GPL570_BR.expr %>%
  as.data.frame() %>%
  mutate(Probe_ID = rownames(.)) %>%
  filter(Probe_ID %in% immune_lncrna) %>% 
  select(Probe_ID,everything())



save(GPL570_BR.expr, GPL570BR_pdata,
     exp_mRNA, exp_lncRNA,
     gpl570_BR_eset.rma,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/2.2gpl570_BR.rdata')


#load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/2.2gpl570_BR.rdata')






