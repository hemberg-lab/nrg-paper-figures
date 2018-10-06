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

logcounts(deng_SCE) <- log2(calculateCPM(deng_SCE, use_size_factors=FALSE) + 1)

deng_SCE <- runPCA(deng_SCE, ntop = 500)

p1 <- ggplot(as.data.frame(reducedDim(deng_SCE)), 
       aes(x = PC1, 
           y = PC2, colour = deng_SCE$cell_type2)) + geom_point() +
  scale_color_tableau("colorblind10") + theme_classic(base_size=12) +
  ggtitle("500 genes") + guides(colour=FALSE)

deng_SCE1 <- deng_SCE[, grepl("8cell", rownames(colData(deng_SCE)))]
deng_SCE1$Protocol <- c(rep("Smart-seq", 28), rep("Smart-seq2", 9))
deng_SCE1 <- runPCA(deng_SCE1, ntop = 500)
p3 <- ggplot(as.data.frame(reducedDim(deng_SCE1)), 
       aes(x = PC1, 
           y = PC2, colour = deng_SCE1$Protocol)) + geom_point() + theme_classic(base_size=12) + guides(colour=guide_legend(title="Protocol"), size = guide_legend(keywidth = 1))

deng_SCE <- runPCA(deng_SCE, ntop = 20000)

p2 <- ggplot(as.data.frame(reducedDim(deng_SCE)), 
       aes(x = PC1, 
           y = PC2, colour = deng_SCE$cell_type2)) + geom_point() +
  scale_color_tableau("colorblind10") + theme_classic(base_size=12) +
  ggtitle("20000 genes") + guides(colour=guide_legend(title="Cell Type"), size = guide_legend(keywidth = 1))

library(cowplot)

plot_grid(p1, p2, ncol = 2, labels = c("a", "b"), label_size = 20, rel_widths = c(1,1.4))
ggsave("pdf/fig2.pdf", w = 9, h = 6)
ggsave("png/fig2.png", w = 9, h = 6)

plot_grid(p3)
ggsave("pdf/fig5.pdf", w = 6, h = 6)
ggsave("png/fig5.png", w = 6, h = 6)
