## 6.1 clusterprofiler

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(cowplot)
library(ReactomePA)
library(tidyverse)
library(enrichplot)
library(DT)
library(pathview)
library(fgsea)


####make gene list from Top Hits limma table
## topTable from make tables script

ENTREZ <- bitr(topTable$gene_name, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Hs.eg.db) ##6% of IDs fail to map
p <- topTable %>% dplyr::select(gene_name, logFC)
p <- p %>% merge(ENTREZ, by.x = "gene_name", by.y = "SYMBOL") %>% distinct()

geneList <- p$logFC
names(geneList) <- as.character(p$ENTREZID)
geneList <- sort(geneList, decreasing = TRUE)
gene <- names(geneList)[abs(geneList) > 2]
geneList <- na.omit(geneList)
head(geneList)

setwd("/scratch/u/rsummey/mcw/data")
write_rds(geneList, "TumorgeneListJune27.rds")
#geneList <- readRDS("TumorgeneListJune27.rds")

#ego <- enrichGO(gene         = names(geneList),
 #                OrgDb         = org.Hs.eg.db,
  #               keyType       = 'SYMBOL',
   #              pAdjustMethod = "BH",
    #             pvalueCutoff  = 0.01,
     #            qvalueCutoff  = 0.05)
#head(ego, 2)

egoCC <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType          = 'ENTREZID',
              ont = 'CC',
              pvalueCutoff = 0.05,
              verbose      = FALSE) 
egoBP <- gseGO(geneList     = geneList,
               OrgDb        = org.Hs.eg.db,
               keyType          = 'ENTREZID',
               ont = 'BP',
               pvalueCutoff = 0.05,
               verbose      = FALSE)
egoMF <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType          = 'ENTREZID',
              ont = 'MF',
              pvalueCutoff = 0.05,
              verbose      = FALSE)
write_rds(egoCC, "clusterGSEAcc.rds")
write_rds(egoBP, "clusterGSEAbpALL_LFCs.rds")
write_rds(egoMF, "clusterGSEAmf.rds")

#egoBP <- read_rds("clusterGSEAbp.rds")
egoCC <- setReadable(egoCC, 'org.Hs.eg.db')
egoBP <- setReadable(egoBP, 'org.Hs.eg.db')
egoMF <- setReadable(egoMF, 'org.Hs.eg.db')

### WRITE TABLES
### SELECT ANYTHING THAT INCLUDES AN ESTROGEN, GnRH, PR, PRL, OR ANDROGEN RECEPTOR

result <- egoCC@result
resultBP <- egoBP@result
resultMF <- egoMF@result

GSEA <- bind_rows(result, resultBP, resultMF)
write.table(GSEA, "GSEAtable_tumorsvsNormaljune27.txt")

datatable(GSEA, escape = FALSE,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'AGCT tumor gene set enrichment analysis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:8), digits=2)

## 14 AR pathways
## 9 PR pathways
## 7 ESR2
## 21 ESR1
## 3 PRLR
## 1 GnRHR1 and 2 (negative)
## 7 LHCGR (all negative)
## 4 FSHR

################ AR table #############
combohormones <- GSEA %>% filter(grepl('/AR/ | /ESR1/ | /ESR2/', core_enrichment))
otherrepro <- GSEA %>% filter(grepl('/PRLR/ | /PGR/ | /GNRH/ | /LHCGR/ | /FSHR/', core_enrichment))


ggplot(combohormones, 
       aes(x=enrichmentScore, y=Description)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="purple", high="orange") +
  theme_bw() +
  labs(title = "GSEA: GO terms including AR, ESR1, ESR2", subtitle ="AGCT Tumors") +
  theme(axis.text.y = element_text(size = 11))

ggplot(otherrepro, 
       aes(x=enrichmentScore, y=Description)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="purple", high="orange") +
  theme_bw() +
  labs(title = "GSEA: GO terms including PRLR, PGR, GNRH, LHCGR, FSHR", subtitle ="AGCT Tumors") +
  theme(axis.text.y = element_text(size = 9))

CCdot <- dotplot(egoCC, showCategory=40) + ggtitle("CC dotplot for GSEA (tumor)")
BPdot <- dotplot(egoBP, showCategory=40) + ggtitle("BP dotplot for GSEA (tumor)")+ 
  theme(axis.text.y = element_text(size = 8))
MFdot <- dotplot(egoMF, showCategory=40) + ggtitle("MF dotplot for GSEA (tumor)")+ 
  theme(axis.text.y = element_text(size = 6))


p1 <- cnetplot(egoCC, categorySize="pvalue", foldChange = geneList) 
p2 <- cnetplot(egoBP, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(egoMF, categorySize="pvalue", foldChange=geneList)

p8 <- cnetplot(egoxBP, categorySize="pvalue", foldChange=list2)
p9 <- cnetplot(egoxMF, categorySize="pvalue", foldChange=list2)


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
#y <- gsePathway(geneList, 
 #               pvalueCutoff = 0.05,
  #              pAdjustMethod = "BH", 
   #             verbose = FALSE)
#dup <- which(duplicated(names(geneList)))
#geneList2 <- geneList[-dup]

#y
#write_rds(y, "reactometumorResults.rds")
#reactomeresult <- readRDS("reactometumorResults.rds")
#viewPathway(
 # "Nuclear Receptor transcription pathway",
  #organism = "human",
  #readable = TRUE,
  #foldChange = geneList2,
  #keyType = "ENTREZID",
  #layout = "kk")
