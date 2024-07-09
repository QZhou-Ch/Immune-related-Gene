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
library(Rsubread)
library(refGenome)
library(GenomicRanges)
library(Rsamtools)

library(tidyverse)

rootdir <- ""


#####1 整理两个免疫治疗的数据集###################

#####1.1Riaz et al., Cell 2017，GSE91061###################
Riaz2017_counts <- read.csv(str_c(rootdir,"/DATABASE/Immunotherapy_trial_studies/data/GSE91061/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz"),
                          row.names=1) %>% 
  as.matrix()

Riaz2017_expr<- limma::voom(Riaz2017_counts)$E %>% 
  t()%>% 
  as.data.frame() %>% 
  mutate(pat = rownames(.)) %>% 
  select(pat,everything())

#设定有效为PR/CR，定义为1
Riaz2017_pdata <- getGEO('GSE91061', destdir = str_c(rootdir,"/DATABASE/GEO/data"),
                         AnnotGPL = F,
                         getGPL = F)[["GSE91061_series_matrix.txt.gz"]]@phenoData@data[,c(1,47)] %>% 
  filter(!(`response:ch1` == 'UNK')) %>% 
  #设定有效为PR/CR，定义为1
  mutate(resp = ifelse(.$`response:ch1` == 'PRCR',1,0),
         pat = str_replace(.$title, '[-]',".")) %>% 
  select(pat,resp)

Riaz2017 <- left_join(Riaz2017_expr,Riaz2017_pdata,by = 'pat') %>% 
  select(pat,resp, everything()) %>% 
  filter(!is.na(resp))

#####1.2 IMvigor210CoreBiologies###################
#加载整理过的IMvigor210数据
load(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/input/IMvigor210CoreBiologies.rdata'))

IMvigor210_expr <- limma::voom(IMvigor210_counts)$E %>% 
  t() %>% 
  as.data.frame() %>% 
  dplyr::mutate(pat = rownames(.)) %>% 
  select(pat,everything()) %>% 
  as.data.frame()

IMvigor210_pdata <- IMvigor210_pData[,c(1,2)] %>% 
  filter(!is.na(binaryResponse)) %>% 
  as.data.frame() %>% 
  mutate(pat = rownames(.),
         resp = ifelse(.$`binaryResponse` == "CR/PR",1,0)) %>% 
  select(pat,resp) %>% 
  as.data.frame()

IMvigor210 <- left_join(IMvigor210_expr,IMvigor210_pdata,by = 'pat') %>% 
  select(pat,resp, everything()) %>% 
  filter(!is.na(resp))


######1.3 合并两个数据集###################
ibc_therapy_ls <- list(Riaz2017,IMvigor210)
names(ibc_therapy_ls) <- c('Riaz2017','IMvigor210')
save(ibc_therapy_ls,file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/1.1.1ibc_therapy_ls.rdata'))
#load the immune gene list
immunegene <- read.delim(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/input/immune_gene.txt'))


######2 logistic回归分析筛选基因###################
logi_result <- ibc_therapy_ls %>% 
  map2(names(.),function(x,y){
    #x = ibc_therapy_ls[[1]]
    #y = names(ibc_therapy_ls)[1]
    gene <- x
    rownames(gene) <- gene$pat
    gene <- gene[-1]
    gene <- cbind(gene$resp, gene[,which(colnames(gene) %in% immunegene$ID)]) %>% 
      rename(resp = `gene$resp`) %>% 
      as.data.frame()
    pvalue <- NULL
    for(i in 2:length(colnames(gene))){
      #i=1117
      a = gene[,i]
      mod <- glm(resp~a,data=gene,family = 'binomial')
      if(is.null(pvalue)){
        pvalue <- cbind(colnames(gene)[i], summary(mod)$coefficients[,"Pr(>|z|)"][2])
      }
      else{
        pvalue <- rbind(pvalue,cbind(colnames(gene)[i], summary(mod)$coefficients[,"Pr(>|z|)"][2]))
      }
    }
    pvalue %<>% 
      as.data.frame() %<>% 
      mutate(set = y)
    colnames(pvalue) <- c('entrez','p_value','set')
    x <- pvalue %>% 
      mutate(entrez = as.character(entrez),
             p_value = as.numeric(as.character(p_value)),
             set = as.character(set))
  })


#type of gene name
keytypes(org.Hs.eg.db)

ibc_therapy_mRNA <- logi_result %>% 
  map(function(x){
    #x = logi_result[[1]]
    entrez <- x %>% 
      #as.data.frame() %>% 
      filter(p_value < 0.05) %>% 
      dplyr::select(entrez)
    x = entrez[,1]
  }) %>% 
  purrr::reduce(intersect) %>% 
  clusterProfiler::bitr(fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = "org.Hs.eg.db") %>% 
  filter(!duplicated(ENTREZID))


write.csv(ibc_therapy_mRNA,file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/2ICI_mRNA.csv'))

ibc_therapy_mRNA <- read.csv(file =
                               str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/2ICI_mRNA.csv'))[,c(2,3)]


save(ibc_therapy_ls, logi_result, ibc_therapy_mRNA, 
     file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/1.1ibc_therapy_mRNA.rdata'))


#load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/1.1ibc_therapy_mRNA.rdata')





