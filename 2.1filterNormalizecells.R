#Filter normalize experimental data

library(edgeR)
library(matrixStats)
library(cowplot)
library(tidyverse)
library(dbplyr)
library(ggplot2)
library(tidyverse)
library(caret)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(data.table)
library(SummarizedExperiment)                                     
library(readr)
library(tibble)
library(dplyr)
library(plyr)
library(gt)
library(rtracklayer)

se <- readRDS("tumorse_ensembl.rds")  ###this is how it was labeled in the readData script
anyNA(se) # check if there are NA values

se_dge <- DGEList(counts = assay(se), genes = rowData(se))
str(se_dge)

dim(se_dge$counts)

prop_expressed <- rowMeans(cpm(se_dge) > 5)
table(prop_expressed)
keep <- prop_expressed > 0.3
table(keep) ## false 25325, true 13870

op <- par(no.readonly = TRUE)
par(mfrow = c(1, 2))
hist(cpm(se_dge, log = TRUE), main = "Unfiltered", xlab = "logCPM")
abline(v = log(1), lty = 2, col = 2)
hist(cpm(se_dge[keep, ], log = TRUE), main = "Filtered", xlab = "logCPM")
abline(v = log(1), lty = 2, col = 2) 
par(op) 
se_dge <- se_dge[keep, , keep.lib.sizes = FALSE]## decrease from 
## 13155 genes
se <- se[keep,]
dim(se_dge)
dim(se) ## same as se_dge

## Now making annotation index file ----
gencode_file <- "/scratch/u/rsummey/Homo_sapiens.GRCh38.107.gtf"
gtf <- import.gff(gencode_file, format = "gtf", genome = "hg38", feature.type = "exon") 
#transcript IDs (ENSEMBL ones) are in gtf if needed.  Also gene_name and transcript_name
#split records by gene to group exons of the same gene                                      
grl <- reduce(split(gtf, elementMetadata(gtf)$gene_name))                                      
gene_lengths <- ldply(grl, function(x) {                                                     
  #sum up the length of individual exons                                                    
  return(c(gene_length = sum(width(x))))                                                
}, .id = "gene_symbol")         

genetype <- unique(elementMetadata(gtf)[, c("gene_name", "gene_id", "gene_biotype")])    
head(genetype)
colnames(genetype)[1] <- "gene_symbol"                                          
gene_lengths <- merge(genetype, gene_lengths, by = "gene_symbol")
head(gene_lengths)

#### if Ensembl version numbers are still in place:
#remove ENSEMBL ID version numbers                                                 
#gene_lengths$ensembl_gene_id <- gsub('\\.[0-9]*', '', gene_lengths$ensembl_gene_id) 
saveRDS(gene_lengths, file = "/scratch/u/rsummey/mcw/data/gene_lengths_prepSingscore_with_ensembl.rds")                  
#gene_lengths <- readRDS("/scratch/u/rsummey/mcw/data/gene_lengths_prepSingscore.rds")

#allocate rownames for ease of indexing ----
rownames(gene_lengths) <- gene_lengths$gene_id
rowData(se)$gene_length <- gene_lengths[rownames(se), "gene_length"]
rowData(se)$gene_biotype <- gene_lengths[rownames(se), "gene_biotype"]

table(rowData(se)$gene_biotype)

#annotate gene lengths for the DGE object
se_dge$genes$length <- gene_lengths[rownames(se_dge), "gene_length"]

gene_annot <- rowData(se)  
gene_annot$gene_name <- row.names(gene_annot)
dim(gene_annot) #write down dimensions 13155x2

unprocessedcounts <- assay(se, withDimnames = FALSE)
unprocesseddf <- as.data.frame(unprocessedcounts)
write.csv(unprocesseddf, "/scratch/u/rsummey/mcw/data/countsraExperiments.csv")

#colnames(unprocesseddf)
dge_tmm = calcNormFactors(se_dge, method = "TMM", refColumn = 39)
dim(dge_tmm)
dim(assay(se))
assay(se, withDimnames = FALSE) <- cpm(dge_tmm, log = TRUE, prior.count = 0.1)
saveRDS(se, file = "/scratch/u/rsummey/mcw/data/ExperimentsSummarizedExpt.rds")                  
counts <- assay(se, withDimnames = FALSE)
head(counts)
## double check rownames and colnames on counts
class(counts)
#countsm <- as.matrix(counts)
countsdf <- as.data.frame(counts)
write.csv(countsdf, "/scratch/u/rsummey/mcw/data/ExperimentCounts.csv")
colSums(countsdf)