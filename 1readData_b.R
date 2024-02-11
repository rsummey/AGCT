#Read in data
#Separate by type

#read in data, counts, filter, normalize ----
library(tidyverse)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(SummarizedExperiment)                                     
library(edgeR)          
library(readr)
library(tibble)
library(dplyr)
library(plyr)
library(gt)
library(singscore)
library(tximport)
library(vsn)

sample_info <- read.csv("/scratch/u/rsummey/mcw/studydesignExperiments24_48.csv")

setwd("/scratch/u/rsummey/mcw/data/kallisto")

#Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
#Tx <- as_tibble(Tx)
#Tx <- dplyr::rename(Tx, target_id = tx_id)
#Tx <- dplyr::select(Tx, "target_id", "gene_name")

#Txi_gene <- tximport(path, 
 #                    type = "kallisto", 
  #                   tx2gene = Tx, 
   #                  txOut = FALSE, #determines whether your data represented at transcript or gene level
    #                 countsFromAbundance = "lengthScaledTPM",
     #                ignoreTxVersion = TRUE)


path <- file.path(sample_info$Sample, "abundance.tsv")
file.exists(path)
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_id")

Txi_gene_symbol <- tximport(path, 
                            type = "kallisto", 
                            tx2gene = Tx, 
                            txOut = FALSE, #determines whether your data represented at transcript or gene level
                            countsFromAbundance = "lengthScaledTPM",
                            ignoreTxVersion = TRUE)

library(edgeR)
library(matrixStats)
library(cowplot)

sampleLabels <- sample_info$Sample
myDGEList <- DGEList(Txi_gene_symbol$counts)
log2.cpm <- cpm(myDGEList, log=TRUE)

cpm <- cpm(myDGEList)
keepers <- rowSums(cpm>5)>=3 #user defined
table(keepers) ##toss 24285 genes, keep 11163
myDGEList.filtered <- myDGEList[keepers,]

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)

myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
counts <- cpm(myDGEList.filtered.norm, log=TRUE) 
write.csv(counts, "experimentCounts_24_48h.csv")
