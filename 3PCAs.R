#PCAs

#PCA ----
library(DT)
library(gt)
library(plotly)
library(tidyverse)

sample_info <- read.table("/scratch/u/rsummey/mcw/studydesignExperiments.txt")
group <- sample_info$group
group <- factor(group)
sample_names <- sample_info$BioSample

### counts from filternormalize

counts <- as.matrix(counts)
counts[is.na(counts)] <- 0

distance <- dist(t(counts), method = "maximum")
clusters <- hclust(distance, method = "average")
plot.new()
par(mfrow = c(1,1))
plot(clusters, labels=sample_names, cex = 0.8)

pca.res <- prcomp(t(counts), scale.=F, retx=T)
ls(pca.res)
summary(pca.res)
pca.res$rotation
pca.res$x
screeplot(pca.res)

pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pc.per
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=sample_names, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot: normalized in cell lines only with C_KGN_1 as reference",
       subtitle = "ADL_KGN_5 outlier",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot)

png(file = "/scratch/u/rsummey/mcw/figures/pcaLinesNewNormalization.png")
pca.plot
dev.off()

pca.plot2 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC3, label=sample_names, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA plot PC3: normalized in cell lines only with C_KGN_1 as reference",
       subtitle = "ADL_KGN_5 and C_070_3 outliers",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot2)

png(file = "/scratch/u/rsummey/mcw/figures/pcaPC3LinesNewNormalization.png")
pca.plot2
dev.off()

#############################remove ADL_KGN_5, C_070_3

sample_info2 <- sample_info[sample_info$BioSample != "ADL_KGN_5" & 
                              sample_info$BioSample != "C_070_3",]
group2 <- sample_info2$group
group2 <- factor(group2)
sample_names2 <- sample_info2$BioSample

counts2 <- counts[,colnames(counts) %in% sample_names2]
counts2[is.na(counts2)] <- 0

distance2 <- dist(t(counts2), method = "maximum")
clusters2 <- hclust(distance2, method = "average")
plot.new()
par(mfrow = c(1,1))
plot(clusters2, labels=sample_names2, cex = 0.8)

pca.res2 <- prcomp(t(counts2), scale.=F, retx=T)
ls(pca.res2)
screeplot(pca.res2)

pc.var2 <- pca.res2$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per2 <- round(pc.var2/sum(pc.var2)*100, 1) 
pc.per2
pca.res.df2 <- as_tibble(pca.res2$x)
pca.plot2 <- ggplot(pca.res.df2) +
  aes(x=PC1, y=PC2, label=sample_names2, color = group2) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot 4",
       subtitle = "5 outliers removed",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot2)

png(file = "/scratch/u/rsummey/mcw/figures/pca4 normalized cells only.png")
pca.plot2
dev.off()

pca.plot3 <- ggplot(pca.res.df2) +
  aes(x=PC1, y=PC3, label=sample_names2, color = group2) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="PCA4 PC3 re normalized",
       subtitle = "5 outliers removed",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot3)

png(file = "/scratch/u/rsummey/mcw/figures/pca4 PC3 normalized cells only.png")
pca.plot3
dev.off()


######################################
## remove C_003_2 and ADL_KGN_1
sample_info2 <- sample_info2[sample_info2$BioSample != "ADL_KGN_1" & 
                              sample_info2$BioSample != "C_003_2",]
group2 <- sample_info2$group
group2 <- factor(group2)
sample_names2 <- sample_info2$BioSample

counts2 <- counts2[,colnames(counts2) %in% sample_names2]
counts2[is.na(counts2)] <- 0

distance2 <- dist(t(counts2), method = "maximum")
clusters2 <- hclust(distance2, method = "average")
plot.new()
par(mfrow = c(1,1))
plot(clusters2, labels=sample_names2, cex = 0.8)

pca.res2 <- prcomp(t(counts2), scale.=F, retx=T)
ls(pca.res2)
screeplot(pca.res2)

pc.var2 <- pca.res2$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per2 <- round(pc.var2/sum(pc.var2)*100, 1) 
pc.per2
pca.res.df2 <- as_tibble(pca.res2$x)
pca.plot2 <- ggplot(pca.res.df2) +
  aes(x=PC1, y=PC2, label=sample_names2, color = group2) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot 4",
       subtitle = "5 outliers removed",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot2)

png(file = "/scratch/u/rsummey/mcw/figures/pca4 normalized cells only.png")
pca.plot2
dev.off()

pca.plot3 <- ggplot(pca.res.df2) +
  aes(x=PC1, y=PC3, label=sample_names2, color = group2) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[3],"%",")")) +
  labs(title="PCA4 PC3 re normalized",
       subtitle = "5 outliers removed",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot3)

png(file = "/scratch/u/rsummey/mcw/figures/pca4 PC3 normalized cells only.png")
pca.plot3
dev.off()


### overwrite counts matrix flat file with outliers removed
countsdf <- as.data.frame(counts2)
write.csv(countsdf, "/scratch/u/rsummey/mcw/data/ExperimentCounts_nooutliers.csv")
counts2 <- read.csv("/scratch/u/rsummey/mcw/data/ExperimentCounts_nooutliers.csv")
write.csv(sample_info2, "/scratch/u/rsummey/mcw/data/studydesignExperiments_nooutliers.csv")
