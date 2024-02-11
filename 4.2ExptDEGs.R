### DEGs anastrozole only

### DEGs by cell line

library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(dplyr)
library(tidyverse)
library(dbplyr)
library(htmltools)

setwd("/scratch/u/rsummey/mcw/data")
sample_info <- read.csv("/scratch/u/rsummey/mcw/data/studydesignExperiments_nooutliers.csv")
counts <- read.csv("/scratch/u/rsummey/mcw/data/ExperimentCounts_nooutliers.csv")
rownames(counts) <- counts$X
counts %<>% dplyr::select(-X)

sample_names <- sample_info$BioSample

category <- factor(sample_info$group, levels = c('c', 'a', 'adl'))
expt <- factor(sample_info$experiment)
line <- factor(sample_info$line)
design <- model.matrix(~0 + category + expt) ##covariate = experiment

dupcor <- duplicateCorrelation(counts, design, block=sample_info$line) ##account for correlation
########################################################################within cell line
dupcor$consensus.correlation
fit <- lmFit(counts, design, block=sample_info$line, correlation=dupcor$consensus.correlation)
contrasts <- makeContrasts(categorya-categoryc, categoryadl-categoryc, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2, trend=TRUE)

a <- topTable(fit3, coef="categorya - categoryc", number = 40000, sort.by = "p")
### no significant DEGs here ^ - minimum adj p value is 0.5
combo <- topTable(fit3, coef="categoryadl - categoryc", number = 40000, sort.by = "p", p.value = 0.05)
### to have expression of all genes for pathway analysis:
combo_not_sig <- topTable(fit3, coef="categoryadl - categoryc", number = 40000, sort.by = "p")

### Tx from read data script

combo$gene_id <- rownames(combo)
combo <- merge(combo, Tx, by = "gene_id", all.y = FALSE) %>% distinct()
combo <- combo %>% dplyr::select(gene_name, gene_id, logFC, adj.P.Val)
datatable(combo,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Differentially expressed genes in cell lines 72 hours after combination treatment',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 25, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:4), digits=2)
write.table(combo, "combination_treatment_DEGs_ensembl.txt")

combo_not_sig$gene_id <- rownames(combo_not_sig)
combo_not_sig <- merge(combo_not_sig, Tx, by = "gene_id", all.y = FALSE) %>% distinct()
combo_not_sig <- combo_not_sig %>% dplyr::select(gene_name, gene_id, logFC, adj.P.Val)
datatable(combo_not_sig,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'All genes in cell lines 72 hours after combination treatment',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 25, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:4), digits=2)
write.table(combo_not_sig, "combination_treatment_all_genes_ensembl.txt")

vplot1 <- ggplot(combo) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", "geneID")) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  labs(title="Top differentially expressed genes in AGCT cell lines and tumors",
       subtitle = "Top DEGs",
       caption=paste0("produced on ", Sys.time())) +
  theme_linedraw()

ggplotly(vplot1)
png("/scratch/u/rsummey/mcw/figures/volcano_DEGs_experiments.png")
vplot1
dev.off()

## No need to do this one since not significant:
#a$gene_id <- rownames(a)
#a <- merge(a, Tx, by = "gene_id", all.y = FALSE) %>% distinct()
#a <- a %>% dplyr::select(gene_name, gene_id, logFC, adj.P.Val)
#datatable(a,
 #         extensions = c('KeyTable', "FixedHeader"),
  #        caption = 'All genes in cell lines 72 hours after anastrozole only treatment',
   #       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 25, lengthMenu = c("10", "25", "50", "100"))) %>%
  #formatRound(columns=c(3:4), digits=2)
#write.table(combo_not_sig, "anastrozole_only_treatment_all_genes_ensembl.txt")

#results <- decideTests(fit3, coef = "categorya - categoryc", method="global")
#results <- as.data.frame(results)
#a<- results[results$`categorya - categoryc` %in% c(-1, 1),]
#a$gene_id <- rownames(a)
#a <- merge(a, Tx, by = "gene_id", all.y = FALSE) %>% distinct()

#write.table(a, "AnastrozoleOnlyGlobalUpDown.txt")

ex <- diffSplice(fit2, geneid=rownames(counts))
## no genes with more than one exon - don't think there are splicing differences here
