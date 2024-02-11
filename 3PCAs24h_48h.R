#PCAs

#PCA ----
library(DT)
library(gt)
library(plotly)
library(tidyverse)
library(dplyr)
library(magrittr)

sample_info <- read.csv("/scratch/u/rsummey/mcw/studydesignExperiments24_48.csv")


head(counts)

group <- sample_info$Condition
group <- factor(group)
sample_names <- sample_info$Sample

colnames(counts) <- sample_names

counts %<>% as.matrix()
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
  labs(title="PCA plot: Experiment timepoints",
       subtitle = "24 and 48 hours",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot)

png(file = "/scratch/u/rsummey/mcw/figures/pca_experimenttimepoints.png")
pca.plot
dev.off()

pca.plot2 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC3, label=sample_names, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC3 (",pc.per[2],"%",")")) +
  labs(title="PCA3 plot: Experiment timepoints",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()
ggplotly(pca.plot2)

png(file = "/scratch/u/rsummey/mcw/figures/pcaPC3_experimenttimepoints.png")
pca.plot2
dev.off()
