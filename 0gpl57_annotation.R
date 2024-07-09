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

rootdir <- ''




# 0.准备工作：构建index ----------------------------------------------------------

if(F){
  #使用R语言构建index
  ref <- str_c(rootdir,
               "/seq_project/reference/ref/ENSEMBL/release_104/Homo_sapiens.GRCh38.dna.primary_assembly.fa")
  buildindex(basename="/seq_project/reference/Rsubread/ENSEMBL/release_104/index",
             reference=ref)
}

# 1.比对生成BAM ---------------------------------------------------------------


#输入芯片探针的fasta文件
reads <- "/seq_project/reference/affymetrix/HG-U133_Plus_2.probe_fasta"

#align
align(index = "/seq_project/reference/Rsubread/ENSEMBL/release_104/index",
      readfile1 = reads,
      output_file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/alignResults.BAM'),
      nthreads = 20, phredOffset=64) 

# 2.将BAM转换成GRanges对象 ------------------------------------------------------

#读取BAM文件
bamFile <- str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/rdata/alignResults.BAM')
quickBamFlagSummary(bamFile)
bam <- scanBam(bamFile)


#将BAM转换成矩阵
tmp <- lapply(bam[[1]], as.character) %>% 
  as.data.frame() %>% 
  #注意flag=4，表示segment unmapped
  filter(flag!=4) %>% 
  mutate(qname = str_split(qname,pattern = ":",simplify = T)[,3])

#将矩阵转换成GRanges对象
my_seq <- with(tmp, GRanges(as.character(rname), #指定染色体号
                            IRanges(as.numeric(pos)-20, as.numeric(pos)+20), #指定范围
                            as.character(strand), #指定strand
                            id = as.character(qname))) #指定id

my_seq <- GRanges(as.character(tmp$rname), #指定染色体号
                  IRanges(as.numeric(tmp$pos)-20, as.numeric(tmp$pos)+20), #指定范围
                  as.character(tmp$strand), #指定strand
                  id = as.character(tmp$qname))#指定id


class(my_seq)


# 3.将参考GTF文件换成GRanges对象 ---------------------------------------------------


gtf <- str_c(rootdir,
             "/seq_project/reference/ref/ENSEMBL/release_104/Homo_sapiens.GRCh38.104.chr.gtf")
#加载GTF文件
ens <- ensemblGenome()
read.gtf(ens,useBasedir = F,gtf)

#将GTF转换成GRanges对象
my_refGTF <- with(getGenePositions(ens), #生成矩阵
                  GRanges(seqid, IRanges(start, end), strand = strand, id = gene_id))


# 4.overlap ---------------------------------------------------------------
#目前有两个GRanges对象，一个是探针my_seq，一个是参考基因my_refGTF，
#需要将探针的位置比对到参考基因的位置上，对应探针号::基因名

overlap <- findOverlaps(my_seq,my_refGTF)

Probe_Ensembl <- cbind(as.data.frame(my_seq[queryHits(overlap)]),
                      as.data.frame(my_refGTF[subjectHits(overlap)]))
Probe_Ensembl <- Probe_Ensembl[,c(6,12)] 
colnames(Probe_Ensembl) <- c("Probe_ID","Ensembl")



# #5.取得唯一的对应即：一个探针只对应一个基因 -------------------------------------------------

Probe_Ensembl <- Probe_Ensembl %>% 
  mutate(temp = str_c(Probe_ID,Ensembl,sep = "_")) %>% 
  filter(!duplicated(temp))

#取得唯一的probe
unique_probe <- setdiff(Probe_Ensembl$Probe_ID,Probe_Ensembl$Probe_ID[duplicated(Probe_Ensembl$Probe_ID)])
#取得唯一的ensemb
unique_ensemb <- setdiff(Probe_Ensembl$Ensembl,Probe_Ensembl$Ensembl[duplicated(Probe_Ensembl$Ensembl)])

#找到一个prob只对应一个基因，但允许一个基因对应到多个prob上
Probe_Ensembl <- Probe_Ensembl %>% 
  filter(Probe_ID %in% unique_probe) %>% 
  select(Probe_ID,Ensembl)

#选择人类数据库
ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
#查看输入输出的数据内容
attributes = listAttributes(ensembl)

#转换
Ensembl <- biomaRt::getBM(
  #指定数据库及物种
  mart = biomaRt::useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl"),
  #指定输入类型
  filters = 'ensembl_gene_id', 
  #指定输出类型
  attributes = c("ensembl_gene_id","gene_biotype"), 
  #设置输入变量
  values = Probe_Ensembl$Ensembl[!duplicated(Probe_Ensembl$Ensembl)],
  useCache = T) %>% 
  rename(Ensembl = ensembl_gene_id) %>% 
  filter()
  

Probe_Ensembl <- left_join(Probe_Ensembl,Ensembl) %>% 
  filter(gene_biotype %in% c("protein_coding","lncRNA")) %>% 
  select(Probe_ID,Ensembl,gene_biotype)






write.csv(Probe_Ensembl,
          file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/probe2ensembl.csv'),
          row.names = F)

Probe_Ensembl <- read.csv(file = str_c(rootdir,'/DATABASE/RProject/ICB_in_BRCA/NEW_ARTICLE/result/probe2ensembl.csv'))










