## 6.1 clusterprofiler

library(clusterProfiler)
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(GOplot)
library(gplots)
library(DT)
library(GSEABase)
library(Biobase)
library(dplyr)
library(GSVA) 

library(gprofiler2)
library(msigdbr)
library(enrichplot)
library(clusterProfiler)
library(limma)
library(cowplot)
library(tidyverse)
library(DOSE)
library(ReactomePA)

hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
hs_gsea %>%
  dplyr::distinct(gs_cat, gs_subcat) %>%
  dplyr::arrange(gs_cat, gs_subcat)

hs_gsea_c5 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C5") %>% # choose your msigdb collection of interest
  dplyr::select(ensembl_gene, entrez_gene)
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(ensembl_gene, entrez_gene)

symbolandentrez <- hs_gsea[,c(5,7)]

### combo from 4.2ExptDEGs script
geneList <- combo$logFC
names <- rownames(combo)
names(geneList) <- names
geneList <- sort(geneList, decreasing = TRUE)
write_rds(geneList, "ComboGeneList.rds")
write.csv(geneList, "ComboGeneList.csv")

geneList2 <- combo$logFC
names <- combo$gene_name
names(geneList2) <- names
geneList2 <- sort(geneList2, decreasing = TRUE)
df <- as.data.frame(names(geneList2))
df$logFC <- geneList2
write.csv(df, "ComboGeneList_gene_symbol.csv")


egoCC <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType          = 'ENSEMBL',
               ont = 'CC',
               pvalueCutoff = 0.05,
               verbose      = FALSE) 
egoBP <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType          = 'ENSEMBL',
               ont = 'BP',
               pvalueCutoff = 0.05,
               verbose      = FALSE)
egoMF <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType          = 'ENSEMBL',
               ont = 'MF',
               pvalueCutoff = 0.05,
               verbose      = FALSE)

#egoBP <- read_rds("clusterGSEAbp.rds")
egoCC <- setReadable(egoCC, 'org.Hs.eg.db')
egoBP <- setReadable(egoBP, 'org.Hs.eg.db')
egoMF <- setReadable(egoMF, 'org.Hs.eg.db')

resultCC <- egoCC@result
resultBP <- egoBP@result
resultMF <- egoMF@result

#datatable(resultBP, escape = FALSE,
 #         extensions = c('KeyTable', "FixedHeader"),
  #        caption = 'AGCT tumor gene set enrichment analysis: BP',
   #       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  #formatRound(columns=c(3:5), digits=2)

#datatable(resultMF, escape = FALSE,
 #         extensions = c('KeyTable', "FixedHeader"),
  #        caption = 'AGCT tumor gene set enrichment analysis: MF',
   #       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  #formatRound(columns=c(3:5), digits=2)

write.csv(result, "comboGSEA_CC.csv")
write.csv(resultBP, "comboGSEA_BP.csv")
write.csv(resultMF, "comboGSEA_MF.csv")


#CCdot <- dotplot(egoCC, showCategory=30) + ggtitle("CC dotplot for GSEA")+ 
 # theme(axis.text.y = element_text(size = 8))
#BPdot <- dotplot(egoBP, showCategory=30) + ggtitle("BP dotplot for GSEA")+ 
 # theme(axis.text.y = element_text(size = 8))
#MFdot <- dotplot(egoMF, showCategory=20) + ggtitle("MF dotplot for GSEA")+ 
 # theme(axis.text.y = element_text(size = 6))

### gene sets significant in each condition
attach(resultCC)
CCtop <- resultCC[order(p.adjust),]
detach(resultCC)
CCtop <- CCtop[1:40,]
CC <- ggplot(CCtop, 
       aes(x=enrichmentScore, y=Description)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="purple", high="orange") +
  theme_bw() +
  labs(title = "GSEA: CC", subtitle ="Combination treatment") +
  theme(axis.text.y = element_text(size = 10))

png(filename = "comboTx_GSEA_CCJun27.png")
CC
dev.off()

attach(resultBP)
BPtop <- resultBP[order(p.adjust),]
detach(resultBP)
BPtop <- BPtop[1:40,]
BP <- ggplot(BPtop, 
       aes(x=enrichmentScore, y=Description)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="purple", high="orange") +
  theme_bw() +
  labs(title = "GSEA: BP", subtitle ="Combination treatment") +
  theme(axis.text.y = element_text(size = 10))
png(filename = "comboTx_GSEA_BPJun27.png")
BP
dev.off()

attach(resultMF)
MFtop <- resultMF[order(p.adjust),]
detach(resultMF)
MFtop <- MFtop[1:40,]
MF <- ggplot(MFtop, 
       aes(x=enrichmentScore, y=Description)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="purple", high="orange") +
  theme_bw() +
  labs(title = "GSEA: MF", subtitle ="Combination treatment") +
  theme(axis.text.y = element_text(size = 10))
png(filename = "comboTx_GSEA_MFJun27.png")
MF
dev.off()


p1 <- cnetplot(egoCC, categorySize="pvalue",foldChange = geneList) 
p2 <- cnetplot(egoBP, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(egoMF, categorySize="pvalue", foldChange=geneList)

x2 <- pairwise_termsim(egoCC)
x3 <- pairwise_termsim(egoBP)
x4 <- pairwise_termsim(egoMF)

p5 <- emapplot(x2, pie = "count") 
p6 <- emapplot(x3, pie = "count")
p7 <- emapplot(x4, pie = "count")

upsetplot(x2)
upsetplot(x3)
upsetplot(x4)


## reactome - pathways down, mostly immune system
############ need entrez for these:
#y <- gsePathway(geneList, 
 #               pvalueCutoff = 0.05,
  #              pAdjustMethod = "BH", 
   #             verbose = FALSE)

#write_rds(y, "reactometumorResults.rds")
#reactomeresult <- readRDS("reactometumorResults.rds")
#viewPathway(
 # "APC/C-mediated degradation of cell cycle proteins",
  #organism = "human",
  #readable = TRUE,
  #foldChange = geneList,
  #keyType = "ENTREZID",
  #layout = "kk")
########### repeat the above with any desired reactome pathway