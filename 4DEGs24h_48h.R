#DEGs

library(limma)
library(edgeR)
library(gt)
library(DT)
library(plotly)
library(dplyr)
library(tidyverse)
library(dbplyr)
library(htmltools)
library(DESeq2)

plotMDS(counts, col = as.numeric(group))
sample_info <- read.csv("/scratch/u/rsummey/mcw/studydesignExperiments24_48.csv")
counts <- read.csv("/scratch/u/rsummey/mcw/data/countsmatrixExperimentTimepoints.csv")
head(counts)
rownames(counts) <- counts$X
counts <- counts %>% dplyr::select(-X)
colnames(counts) <- sample_info$Sample

category <- factor(sample_info$Condition, levels = c('Combo', 'Anastrozole', 'Control'))
time <- factor(sample_info$Time)
experiment <- factor(sample_info$Experiment)
design <- model.matrix(~0 + category + time + experiment)
colnames(design) <- c("Combo", "Anastrozole", "Control", "Time48", "Experiment2")

counts[counts<0] <- 0

v.counts <- voom(counts, design, plot = TRUE)
fit <- lmFit(v.counts, design)

contrast.matrix1 <- makeContrasts(ComboTx = Combo - Control,
                                 levels = colnames(coef(fit)))
contrast.matrix2 <- makeContrasts(AIonly = Anastrozole - Control,
                                levels = colnames(coef(fit)))

contrastfit1 <- contrasts.fit(fit, contrast.matrix1)
contrastfit2 <- contrasts.fit(fit, contrast.matrix2)

myebFit <- eBayes(fit)
ebContrast1 <- eBayes(contrastfit1)
ebContrast2 <- eBayes(contrastfit2)

myTopHits <- topTable(myebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC",
                      p.value = 0.05)
length(which(myTopHits$adj.P.Val < 0.05))
#no of DEGs - 11484

contrastTopHits1 <- topTable(ebContrast1, adjust ="BH", coef=1, number=40000, sort.by = "logFC", 
                             p.value=0.05)
length(which(contrastTopHits1$P.Val < 0.05))
#no of DEGs - 1008

contrastTopHits2 <- topTable(ebContrast2, adjust ="BH", coef=1, sort.by = "logFC",
                             number = "inf",
                             p.value = 0.05)
length(which(contrastTopHits2$adj.P.Val < 0.05))
#no of DEGs - 0


#######################timepoints individually
#######################
sample_info24 <- sample_info[sample_info$Time == 24,]
sample_info48 <- sample_info[sample_info$Time == 48,]
counts24 <- counts[,colnames(counts) %in% sample_info24$Sample]
counts48 <- counts[,colnames(counts) %in% sample_info48$Sample]

category24 <- factor(sample_info24$Condition, levels = c('Combo', 'Anastrozole', 'Control'))
category48 <- factor(sample_info48$Condition, levels = c('Combo', 'Anastrozole', 'Control'))

experiment24 <- factor(sample_info24$Experiment)
experiment48 <- factor(sample_info48$Experiment)

design24 <- model.matrix(~0 + category24 + experiment24)
design48 <- model.matrix(~0 + category48 + experiment48)
colnames(design24) <- c("Combo", "Anastrozole", "Control", "experiment242")
colnames(design48) <- c("Combo", "Anastrozole", "Control", "experiment482")

v.counts24 <- voom(counts24, design24, plot = TRUE)
v.counts48 <- voom(counts48, design48, plot = TRUE)

fit24 <- lmFit(v.counts24, design24)
fit48 <- lmFit(v.counts48, design48)


contrast.matrix124 <- makeContrasts(ComboTx = Combo - Control,
                                  levels = colnames(coef(fit24)))
contrast.matrix148 <- makeContrasts(ComboTx = Combo - Control,
                                    levels = colnames(coef(fit48)))
contrast.matrix224 <- makeContrasts(AIonly = Anastrozole - Control,
                                  levels = colnames(coef(fit24)))
contrast.matrix248 <- makeContrasts(AIonly = Anastrozole - Control,
                                    levels = colnames(coef(fit48)))

contrastfit1_24 <- contrasts.fit(fit24, contrast.matrix124)
contrastfit1_48 <- contrasts.fit(fit48, contrast.matrix148)

contrastfit2_24 <- contrasts.fit(fit24, contrast.matrix224)
contrastfit2_48 <- contrasts.fit(fit48, contrast.matrix248)

myebFit24 <- eBayes(fit24)
myebFit48 <- eBayes(fit48)

ebContrast1_24 <- eBayes(contrastfit1_24)
ebContrast1_48 <- eBayes(contrastfit1_48)

ebContrast2_24 <- eBayes(contrastfit2_24)
ebContrast2_48 <- eBayes(contrastfit2_48)


myTopHits24 <- topTable(myebFit24, adjust ="BH", coef=1, number=40000, sort.by="logFC",
                      p.value = 0.05)
length(which(myTopHits24$adj.P.Val < 0.05))
#no of DEGs - 11480
myTopHits48 <- topTable(myebFit48, adjust ="BH", coef=1, number=40000, sort.by="logFC",
                        p.value = 0.05)
length(which(myTopHits48$adj.P.Val < 0.05))

contrastTopHits1_24 <- topTable(ebContrast1_24, adjust ="BH", coef=1, number=40000, sort.by = "logFC", 
                             p.value=0.05)
length(which(contrastTopHits1_24$P.Val < 0.05))
contrastTopHits1_48 <- topTable(ebContrast1_48, adjust ="BH", coef=1, number=40000, sort.by = "logFC", 
                                p.value=0.05)
length(which(contrastTopHits1_48$P.Val < 0.05))
#no of DEGs - 1

contrastTopHits2_24 <- topTable(ebContrast2_24, adjust ="BH", coef=1, sort.by = "logFC", number = "inf")
length(which(contrastTopHits2_24$adj.P.Val < 0.05))
contrastTopHits2_48 <- topTable(ebContrast2_48, adjust ="BH", coef=1, sort.by = "logFC", number = "inf")
length(which(contrastTopHits2_48$adj.P.Val < 0.05))
#no of DEGs - 0

