### tumor z scores

library(magrittr)
library(dplyr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(DT)

setwd("/scratch/u/rsummey/mcw/data")
counts <- read.table("tumorandnormalovarycounts.txt")
head(counts[,1:10])
rownames(counts) <- counts$Row.names
counts %<>% dplyr::select(-Row.names)
z <- scale(counts)
head(z[,1:8])
z <- z[,1:9]
z %<>% as.data.frame()
z$ensemblID <- rownames(z)

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::select(Tx, "gene_name", "gene_id")

z %<>% as.data.frame()
class(z)
z <- merge(Tx, z, by.x = "gene_id", by.y = "ensemblID", all = FALSE) %>% distinct()
z <- z %>% mutate_at(c('T2012_020', 'T2018_002', 'T2011_070_2', 'T2012_060_2',
                       'T2015_071_2', 'T2018_004_2', 'T2013_192',
                       'T2022_032', 'T2013_005_2'), as.numeric)

z <- z %>% rowwise() %>% 
  mutate(avgZ = mean(c_across(c('T2012_020', 'T2018_002', 'T2011_070_2', 'T2012_060_2',
                                         'T2015_071_2', 'T2018_004_2', 'T2013_192',
                                         'T2022_032', 'T2013_005_2')), na.rm = TRUE))
prop_expressed <- rowMeans(z > 1.96)
table(prop_expressed)
keep <- prop_expressed > 0.3
sig <- z[keep,]
sig %<>% na.omit()
write.table(sig, "zscore_sig_in_3_jun28.txt")

datatable(z, escape = FALSE,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'AGCT Gene z-scores',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:12), digits=3)

degs <- read.table("degsTumorLogfc.txt")
head(degs)
inter <- intersect(sig$gene_name, degs$gene_name)

deg_z <- degs[degs$gene_name %in% c(inter),]
write.table(deg_z, "deg_zscoreintersect.txt")
datatable(deg_z, escape = FALSE,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'AGCT intersection of DEGs and z-scores',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:4), digits=4)
