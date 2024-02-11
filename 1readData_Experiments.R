#Read in data from experimental replicates

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

sample_info <- read.table("/scratch/u/rsummey/mcw/studydesignExperiments.txt")


setwd("/scratch/u/rsummey/mcw/data/kallisto")

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "gene_name", "entrez_id")) ## pull transcript maps from Ensembl
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_id", "gene_name")

path <- file.path(sample_info$BioSample, "abundance.tsv") ## set path to files
file.exists(path)


kallisto2se <- function(sample_info = sample_info, level=level, ###make a custom function to read kallisto counts into a summarized experiment
                        kallistodir = kallistodir, 
                        tx2gene = tx2gene,
                        countsFromAbundance = "lengthScaledTPM") {
  
  sample_meta <- sample_info[!duplicated(sample_info$BioSample), ]
  files <- file.path(kallistodir)
  names(files) <- paste0(sample_meta$BioSample)
  coldata <- data.frame(row.names = sample_meta$BioSample)
  coldata <- cbind(coldata, sample_meta[, !names(sample_meta) == "BioSample"])
  
  if(level == "gene") {
    exp <- tximport::tximport(
      files, type = "kallisto", tx2gene = tx2gene,
      ignoreTxVersion = TRUE,
      countsFromAbundance = countsFromAbundance
    )
    final <- SummarizedExperiment::SummarizedExperiment(
      assays = list(gene_TPM = exp$abundance, gene_counts = exp$counts),
      colData = coldata
    )
  } else if(level == "transcript") {
    exp <- tximport::tximport(
      files, type = "kallisto", txOut = TRUE,
      ignoreTxVersion = TRUE,
      countsFromAbundance = countsFromAbundance
    )
    final <- SummarizedExperiment::SummarizedExperiment(
      assays = list(tx_TPM = exp$abundance, tx_counts = exp$counts),
      colData = coldata
    )
  } else if(level == "both") {
    exp_tx <- tximport::tximport(
      files, type = "kallisto", txOut = TRUE,
      ignoreTxVersion = TRUE,
      countsFromAbundance = countsFromAbundance
    )
    exp_gene <- tximport::summarizeToGene(exp, tx2gene)
    se_gene <- SummarizedExperiment::SummarizedExperiment(
      assays = list(gene_TPM = exp_gene$abundance, 
                    gene_counts = exp_gene$counts),
      colData = coldata
    )
    exp_tx <- tximport::summarizeToGene(exp_tx, tx2gene)
    se_tx <- SummarizedExperiment::SummarizedExperiment(
      assays = list(tx_TPM = exp_tx$abundance, 
                    tx_counts = exp_tx$counts),
      colData = coldata
    )
    final <- list(gene = se_gene, transcript = se_tx)
  } else {
    stop("Invalid parameter for the 'level' argument.")
  }
  return(final)
}

##read in data from paths
se <- kallisto2se(sample_info = sample_info, level = "gene", kallistodir = path, 
                       tx2gene = Tx)
str(se)
setwd("/scratch/u/rsummey/mcw/data")
write_rds(se, "tumorse_ensembl.rds")

se <- readRDS("tumorse_ensembl.rds")
mcw <- assay(se)
head(mcw)
