library(SingleCellExperiment)
library(TSCAN)
library(destiny)
library(scater)
library(ggplot2)
library(ggthemes)
library(ggbeeswarm)
library(corrplot)
set.seed(1)

deng_SCE <- readRDS("deng-reads.rds")
deng_SCE$cell_type2 <- as.character(deng_SCE$cell_type2)
deng_SCE$cell_type2[deng_SCE$cell_type2 == "zy"] <- "Zygote"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "early2cell"] <- "2 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "mid2cell"] <- "2 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "late2cell"] <- "2 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "4cell"] <- "4 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "8cell"] <- "8 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "16cell"] <- "16 cells"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "earlyblast"] <- "Blastocyst"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "midblast"] <- "Blastocyst"
deng_SCE$cell_type2[deng_SCE$cell_type2 == "lateblast"] <- "Blastocyst"
deng_SCE$cell_type2 <- factor(
    deng_SCE$cell_type2,
    levels = c("Zygote", "2 cells", "4 cells", "8 cells", "16 cells", "Blastocyst")
)
cellLabels <- deng_SCE$cell_type2

logcounts(deng_SCE) <- log2(calculateCPM(deng_SCE, use_size_factors=FALSE) + 1)

deng <- counts(deng_SCE)
colnames(deng) <- cellLabels

procdeng <- TSCAN::preprocess(deng)
colnames(procdeng) <- 1:ncol(deng)
dengclust <- TSCAN::exprmclust(procdeng, clusternum = 10)
TSCAN::plotmclust(dengclust)
dengorderTSCAN <- TSCAN::TSCANorder(dengclust, orderonly = FALSE)
pseudotime_order_tscan <- as.character(dengorderTSCAN$sample_name)
deng_SCE$pseudotime_order_tscan <- NA
deng_SCE$pseudotime_order_tscan[as.numeric(dengorderTSCAN$sample_name)] <- 
    -dengorderTSCAN$Pseudotime

p1 <- ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_order_tscan, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau("colorblind10") + theme_classic(base_size = 12) +
    xlab("Pseudotime") + ylab("Timepoint") +
  ggtitle("TSCAN") + guides(colour=FALSE)

deng <- logcounts(deng_SCE)
colnames(deng) <- cellLabels
dm <- DiffusionMap(t(deng))

tmp <- data.frame(DC1 = eigenvectors(dm)[,1],
                  DC2 = eigenvectors(dm)[,2],
                  Timepoint = deng_SCE$cell_type2)

deng_SCE$pseudotime_diffusionmap <- rank(eigenvectors(dm)[,1])
p2 <- ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = pseudotime_diffusionmap, 
           y = cell_type2, colour = cell_type2)) +
    geom_quasirandom(groupOnX = FALSE) +
    scale_color_tableau("colorblind10") + theme_classic(base_size = 12) +
    xlab("Pseudotime") +
    ylab("Timepoint") +
  ggtitle("Diffusion Map") + guides(colour=FALSE)

library(cowplot)
plot_grid(p1, p2, ncol = 2, labels = c("a", "b"), label_size = 20)
ggsave("pdf/fig4.pdf", w = 9, h = 6)
ggsave("png/fig4.png", w = 9, h = 6)

