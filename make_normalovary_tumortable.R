#### GET NORMAL OVARY DATA

library(data.table)
library(dplyr)
library(stringr)
library(biomaRt)
library(edgeR)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DT)
library(rDGIdb)
library(EnsDb.Hsapiens.v86)
library(tximport)
library(magrittr)

setwd("/scratch/u/rsummey/mcw/data")

Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, target_id, gene_id, gene_name)

gtex <- fread("/scratch/g/ehopp/work/rnaseq/panca_raw/gtex_Kallisto_est_counts.txt")
gtexmeta <- fread("/scratch/g/ehopp/work/rnaseq/panca_raw/GTEX_phenotype.txt")

gtexmeta <- gtexmeta[gtexmeta$`_primary_site` == "Ovary",]
gtexsamples <- intersect(gtexmeta$Sample, colnames(gtex))
gtex <- gtex %>% dplyr::select(c(sample, gtexsamples))
geneIDs <- gtex$sample
geneIDs %<>% str_remove( "\\..*")
gtex$sample <- geneIDs

gtex <- merge(Tx, gtex, by.x = "target_id", by.y = "sample", 
                    all = FALSE) %>% 
  distinct()  ##gtex comes with target id, labeled in matrix as "sample"
gtex %<>% dplyr::select(-c(gene_name, target_id))

gtex %<>% as.data.table()
gtex2 <- gtex[ , lapply(.SD, sum), by = gene_id]   # Aggregate data
#remove(gtex)

gtex2 %<>% as.data.frame()
rownames(gtex2) <- gtex2$gene_id
gtex2 %<>% dplyr::select(-gene_id)
gtex2 %<>% as.matrix()

#write.table(gtex2, "gtex_ensembl_geneid_reverseLog.txt")

#gtex <- read.table("gtex_ensembl_geneid.txt", row.names = 1)  if you need to read saved data back in

### read MCW tumor samples -----
sample_info <- read.csv("/scratch/u/rsummey/mcw/tumoronlyNoOutliers.csv")
names <- sample_info$BioSample
setwd("/scratch/u/rsummey/mcw/data/kallisto")
path <- file.path(sample_info$BioSample, "abundance.tsv") # set file paths to your mapped datafile.exists(path)
file.exists(path)
Txi_gene <- tximport(path, 
                    type = "kallisto", 
                   tx2gene = Tx, 
                  txOut = FALSE, #determines whether your data represented at transcript or gene level
                 countsFromAbundance = "lengthScaledTPM",
                ignoreTxVersion = TRUE)
agct <- DGEList(Txi_gene$counts)
agct <- agct$counts
colnames(agct) <- names
head(agct)
agct <- agct + 1
agct <- log2(agct)

m <- merge(agct, gtex2, by = 0, all = FALSE) %>% distinct() ###all=FALSE to only include genes analyzed in both
head(m[,1:15])
rownames(m) <- m$Row.names
m %<>% dplyr::select(-Row.names)
setwd("/scratch/u/rsummey/mcw/data")
#write.table(m, "tumorandnormalovaryNon_Normcounts.txt")
#m1 <- read.table("tumorandnormalovarycounts.txt", row.names = 1)

#rownames(m) <- m$Row.names
#m %<>% dplyr::select(-Row.names)
m[is.na(m)] <- 0
m[m<0] <- 0

group <- factor(c(rep("tumor", times = 9), rep("normal", times = 88)))
study_info <- data.frame(colnames(m), group)
design <- model.matrix(~0 + group)

dge <- DGEList(counts=m)
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge)
m <- cpm(dge, log=TRUE, prior.count=3)

#fwrite(m, "normalized_tumor_ovary.csv")
#m <- fread("normalized_tumor_ovary.csv")
#m %<>% as.data.frame()
#rownames(m) <- m$geneID
#m %<>% dplyr::select(-geneID)


plotMDS(m, col = as.numeric(group)) ### exclude GTEX-S4UY

### IF YOU WANT TO USE AGE:
#URL <- "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
##### match to age - postmenopausal and premenopausal groups
#download.file(URL, "gtex_metadata.txt")
#meta <- fread("gtex_metadata.txt")
#head(meta)
#meta <- meta[meta$SUBJID %in% names,]

##### remove subjects without age info (52 left)
##### add age info to design
#meta <- as.data.frame(meta)
#meta <- meta %>% dplyr::select(SUBJID, AGE)
#design
#design <- merge(design, meta, by.x = "names", by.y = "SUBJID", all = TRUE)
## these are sorted by age category so can remove top ones (normal NAs)
#design <- design[36:97,]
#design$AGE[is.na(design$AGE)] <- "AGCT"
#design$AGE[design$AGE %in% c("20-29", "30-39", "40-49")] <- "PRE"
#design$AGE[design$AGE %in% c("50-59", "60-69")] <- "POST"

################ PCAS##################################
colnames(design) <- c("normal", "tumor")
v.counts <- voom(dge, design, plot = TRUE)
fit <- lmFit(v.counts, design)
myebFit <- eBayes(fit)
topHits <- topTable(myebFit, adjust ="BH", coef=1, sort.by = "logFC", number = "inf",
                             p.value = 0.05)

Tx %<>% dplyr::select(-target_id)
topHits <- merge(topHits, Tx, by.x = 0, by.y = "gene_id")
topTable <- topHits %>% dplyr::select(gene_name, logFC, AveExpr, adj.P.Val)
topTable %<>% distinct()
length(unique(topTable$gene_name))

write.table(topTable, "degsTumorLogfc.txt")

datatable(topTable,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Differentially expressed genes in AGCT compared with normal ovary',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 50, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:4), digits=2)

topTable %<>% as_tibble()

vplot1 <- ggplot(topTable) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", gene_name)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 10, ymin = -log10(0.01), ymax = 210, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -2, xmax = -1, ymin = -log10(0.01), ymax = 210, alpha=.2, fill="#2C467A") +
  labs(title="Top differentially expressed genes in AGCT compared with normal ovary",
       subtitle = "Top DEGs",
       caption=paste0("produced on ", Sys.time())) +
  theme_linedraw()

png(filename = "finalVolcanoTumorNormal.png")
vplot1
dev.off()

results <- decideTests(myebFit, method="global", adjust.method="BH", 
                       p.value=0.05, lfc=0) #this is just a yes/no disease
as.table(summary(results)) ############# up in all

list <- topTable$gene_name
write_rds(list, "tumorDEGnamesList.rds")
write.csv(list, "tumorDEGlist")

