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
library(magrittr)
library(furrr)





library(tidyverse)

#####1 做两组患者的差异表达基因################################################
##首先整理患者的pdata
load('/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/3.1GPL570_train.rdata')
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/4.1TCGA_Score.rdata')

tcga_clust <- tcga_pscore[,c(2,4)]

##然后整理患者的counts
load(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/2.1TCGA_counts.rdata')
tcga_rnas <- tcga_BR[,c(2,68:ncol(tcga_BR))] %>% 
  as.data.frame()
rownames(tcga_rnas) <- tcga_rnas$sample
tcga_rnas <- tcga_rnas[which(rownames(tcga_rnas) %in% tcga_clust$sample),] %>% 
  select(-'sample') %>% 
  t()

#可以跳过分析，直接读取结果
Allgene <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/TCGA_Alldegene.txt", row.names=1)



if(F){
  ##执行差异表达分析
  register(MulticoreParam(6))
  dds <- DESeqDataSetFromMatrix(countData = tcga_rnas,
                                colData = tcga_clust,
                                design = ~ Group)
  
  #设定Low_Score组为参照,过滤低表达数据
  dds$Group <- relevel(dds$Group, ref = "Low_Score")
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]
  
  #得到差异基因
  Allgene <- as.data.frame(results(DESeq(dds, parallel = TRUE)), 
                           contrast = c('Group','High_Score','Low_Score'))
  #检查是否有NA值
  table(is.na(Allgene$padj))
  
  write.table(Allgene, "/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/TCGA_Alldegene.txt",sep = '\t',quote = F)
  #Allgene <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/TCGA_Alldegene.txt", row.names=1)
}




Gene_fc <- Allgene %>% 
  filter(padj < 0.1) %>% 
  mutate(gene = rownames(.)) %>% 
  .[order(.$log2FoldChange,decreasing = T),c(2,7)] %>% 
  as.data.frame()

Gene_list <- Gene_fc$log2FoldChange
names(Gene_list) <- Gene_fc$gene

#保存差异表达基因
Diffgene <- Allgene[(Allgene$padj < 0.01 & (abs(Allgene$log2FoldChange) > 1)),]
length(Diffgene$log2FoldChange[which(Diffgene$log2FoldChange > 0)])
length(Diffgene$log2FoldChange[which(Diffgene$log2FoldChange < 0)])
write.table(Diffgene, "/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/TCGA_Diffgene.txt",sep = '\t',quote = F)
#Diffgene <- read.delim("/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/table/TCGA_Diffgene.txt", row.names=1)




#####2 绘制差异表达基因的散点图################################################
TCGA_DE_DOT <- Allgene %>% 
  select(c(2,5,6)) %>% 
  mutate(Group = "non-significant",logp = -log10(padj),ENSEMBL = rownames(.)) %>% 
  as.tibble() %>% 
  left_join(Gene_coef[c(1,2)],by = 'ENSEMBL')

TCGA_DE_DOT[(which((TCGA_DE_DOT$padj < 0.05) & (TCGA_DE_DOT$log2FoldChange > 1))),4] = "up-regulated"
TCGA_DE_DOT[(which((TCGA_DE_DOT$padj < 0.05) & (TCGA_DE_DOT$log2FoldChange < -1))),4] = "down-regulated"
table(TCGA_DE_DOT$Group)


F5A_TCGA_DE_DOT_Plot <- ggscatter(TCGA_DE_DOT,x = "log2FoldChange",y = 'logp',#指定数据集，x轴，y轴
                   color = 'Group',#颜色按照差异低表达，无差异表达，差异高表达分组
                   palette = c("#619cff","#d7d7d8","#f8766d"),#指定分组对应的颜色
                   size = 0.5,
                   #shape = 'Type',#形状按照mRNA，lncRNA，miRNA
                   #alpha = "Expression",#透明度按照表达量展示
                   #size=DotPlot$Expression,#指定点的大小
                   #label = DotPlot$genename,
                   #label.select = DotPlot$genename[!is.na(DotPlot$genename)],#展示前10个差异基因的基因名
                   font.label = 8,
                   repel = T)+
  theme_classic2()+
  geom_hline(yintercept=1.30,linetype ="dashed")+
  geom_vline(xintercept=c(-1,1),linetype ="dashed")+
  labs(title="Different Expression between High_Score and Low_Score", x="Log2 FoldChange", y="-Log10(Adj P-value)")+
  # #增加一个图层
  # ggrepel::geom_label_repel(aes(label = genename),data = TCGA_DE_DOT[-which(is.na(TCGA_DE_DOT$genename)),],
  #                           max.overlaps = getOption("ggrepel.max.overlaps", default = 100),)+
  theme(plot.title = element_text(hjust = 0.5))
F5A_TCGA_DE_DOT_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_A.pdf',width=7,height=7)
F5A_TCGA_DE_DOT_Plot
dev.off()

#####3 根据差异表达基因绘制GO、KEGG富集分析图################################################

##1首先做GSEA分析
immune_cell_matrix <- read.csv(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/result/immune_cell_signature.csv')[,-1] %>% 
  filter(!(Type %in% c('Neutrophils','Gamma_delta_T_cell')))

#GSEA分析并绘图
immune_cell <- GSEA(Gene_list,TERM2GENE = immune_cell_matrix[,c(1,3)])
F5B_immune_cell_GSEA_Plot <- enrichplot::gseaplot2(immune_cell,1:4)
F5B_immune_cell_GSEA_Plot

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_B.pdf',width=7,height=7)
F5B_immune_cell_GSEA_Plot
dev.off()


#GO富集分析
ego <- enrichGO(gene = rownames(Diffgene),
                OrgDb = org.Hs.eg.db,
                keyType = 'ENSEMBL',
                ont = "all",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

ego@result$Description[which(rownames(ego@result)=='GO:0002460')] <- "adaptive immune response based on somatic..."
ego@result <- ego@result[order(ego@result$p.adjust,decreasing = F),]


F5C_enrichGO_Plot <- dotplot(ego, color = "p.adjust", title = "GO Enrichment",orderBy = "x",showCategory= 25)+ 
  facet_grid(ONTOLOGY~.,scale="free",space="free")

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_C.pdf',width=9,height=9)
F5C_enrichGO_Plot
dev.off()



#KEGG富集分析

gene_df <- read.delim("/Volumes/DATA/Annotation/MSigDB/7.2/msigdb_v7.2_chip_files/Human_ENSEMBL_Gene_ID_MSigDB.v7.2.chip")[,-3] %>% 
  left_join(read.delim("/Volumes/DATA/Annotation/MSigDB/7.2/msigdb_v7.2_chip_files/Human_NCBI_Entrez_Gene_ID_MSigDB.v7.2.chip")[,-3],
            by = 'Gene.Symbol') %>% 
  rename('ENSEMBL' = "Probe.Set.ID.x",
         'SYMBOL' = "Gene.Symbol",
         'ENTREZID' = "Probe.Set.ID.y") %>% 
  filter(ENSEMBL %in% rownames(Diffgene))

ekegg <- enrichKEGG(gene = gene_df$ENTREZID,
                    organism = "hsa",
                    keyType = "kegg",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ekegg@result <- ekegg@result[order(ekegg@result$p.adjust,decreasing = F),]
F5D_enrichKEGG_Plot <- dotplot(ekegg,color = "p.adjust", title = "KEGG Enrichment", orderBy = "x",showCategory= 25)

pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_D.pdf',width=9,height=9)
F5D_enrichKEGG_Plot
dev.off()

#####4 基于FPKM，使用ssGSEA完成MSigDB H分析################################################

tcga_FPKM <- tcga_rnaseq_FPKM %>% 
  filter(sample %in% tcga_clust$sample) %>% 
  as.data.frame() 
rownames(tcga_FPKM) <- tcga_FPKM$sample


tcga_FPKM <- tcga_FPKM[,-1] %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Probe.Set.ID = rownames(.)) %>% 
  left_join(read.delim("/Volumes/DATA/Annotation/MSigDB/7.2/msigdb_v7.2_chip_files/Human_ENSEMBL_Gene_ID_MSigDB.v7.2.chip")[,-3],
            by = 'Probe.Set.ID') %>%
  select(-'Probe.Set.ID') %>% 
  filter(!is.na(Gene.Symbol)) %>% 
  select("Gene.Symbol",everything()) %>% 
  dplyr::mutate(median = apply(.[,c(2:ncol(.))],1,median)) %>%
  .[order(.$Gene.Symbol,.$median,decreasing = T),] %>%
  filter(!duplicated(Gene.Symbol)) %>%
  select(-'median') %>% 
  as.data.frame()


rownames(tcga_FPKM) <- tcga_FPKM$Gene.Symbol
tcga_FPKM <- tcga_FPKM[,-1]


ssGSEA_H <- GSVA::gsva(as.matrix(tcga_FPKM), 
                       getGmt("/Volumes/DATA/Annotation/MSigDB/7.2/h.all.v7.2.symbols.gmt", 
                              collectionType=BroadCollection(category="h"), sep="\t"), 
                       method = "ssgsea", ssgsea.norm = F, parallel.sz = 6) %>% 
  as.data.frame() %>% 
  apply(., 1, function(x){
    a = as.numeric(x)
    a = scale(x,center = T,scale = T)
    names(a) = names(x)
    return(a)
  }) %>% 
  as.data.frame() %>% 
  mutate(sample = rownames(.)) %>% 
  left_join(tcga_pscore[,c(2:4)],by = 'sample') %>% 
  as.data.frame() %>% 
  select('sample','Group','score',everything())

ssGSEA_H <- cbind(ssGSEA_H$sample,ssGSEA_H$Group,apply(ssGSEA_H[,c(3:ncol(ssGSEA_H))],c(1,2),as.numeric)) %>% 
  as.data.frame()
colnames(ssGSEA_H)[1:2] <- c('sample','Group')

ssGSEA_H <- ssGSEA_H[order(ssGSEA_H$score ,decreasing = T),]




MSigDB_Heat <- ssGSEA_H[,-c(2,3)]

rownames(MSigDB_Heat) <- MSigDB_Heat[,1]
MSigDB_Heat <- t(apply(as.matrix(MSigDB_Heat[,-1]),c(1,2),as.numeric))
rownames(MSigDB_Heat) <- str_split(rownames(MSigDB_Heat),pattern = 'HALLMARK_',simplify = T)[,2]



ANNO_Heat <- ssGSEA_H[,c(1:3)]
rownames(ANNO_Heat) <- ANNO_Heat[,1]


immune_cell <- GSEA(Gene_list,TERM2GENE = immune_cell_matrix[,c(1,3)])
F5B_immune_cell_GSEA_Plot <- enrichplot::gseaplot2(immune_cell,1:4)
F5B_immune_cell_GSEA_Plot






F5E_heatmap_Plot <- Heatmap(MSigDB_Heat,
                        #heatmap_width = unit(10, "cm"), heatmap_height = unit(10, "cm"),
                        width = unit(6, "cm"), height = unit(15, "cm"),
                        row_gap = unit(2, "mm"),
                        heatmap_legend_param = list(title = "NES"),
                        col = circlize::colorRamp2(c(-3, 0, 3), c("#004977", "white", "#d03027")),
                        cluster_rows = T, cluster_columns = F,
                        row_names_gp = gpar(fontsize = 7),
                        show_row_dend = F,
                        row_km = 2, row_title = NULL,
                        show_column_names = F, row_names_side = "left",row_dend_side = "right",
                        bottom_annotation = HeatmapAnnotation(Score = as.numeric(ANNO_Heat$score),
                                                              col = list(Score = circlize::colorRamp2(c(min(as.numeric(ANNO_Heat$score)),
                                                                                                        max(as.numeric(ANNO_Heat$score))),
                                                                                                      c('#ffffff','#338d11')))),
                        top_annotation = HeatmapAnnotation(Group = factor(ANNO_Heat$Group,levels = c('High_Score','Low_Score')),
                                                           col = list(Group = c('High_Score'="#d03027",'Low_Score'="#004977"))))


pdf(file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/figure/F5_E.pdf',width=8,height=8)
F5E_heatmap_Plot
dev.off()


save(F5A_TCGA_DE_DOT_Plot,TCGA_DE_DOT,
     F5B_immune_cell_GSEA_Plot,immune_cell,
     F5C_enrichGO_Plot,ego,
     F5D_enrichKEGG_Plot,ekegg,
     F5E_heatmap_Plot,MSigDB_Heat,ANNO_Heat,
     file = '/Volumes/DATA/RProject/ICB_in_BRCA/MicroArray/rdata/F5A_E.rdata')







