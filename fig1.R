library(SingleCellExperiment)
library(scater)
library(cowplot)
library(statmod)
library(ggthemes)
library(pheatmap)
source("utils.R")

base_font_size = 10

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

x <- calculateQCMetrics(deng_SCE)
logcounts(x) <- log2(calculateCPM(x, use_size_factors=FALSE) + 1)
assay(x, "logcounts_raw") <- log2(counts(x) + 1)
assay(x, "norm") <- calculateCPM(x, use_size_factors=FALSE)


p1 <- ggplot(as.data.frame(colData(deng_SCE)), 
       aes(x = total_features)) +
  geom_histogram() +
  theme_classic(base_size = base_font_size) +
  geom_vline(xintercept = 4500, colour = "red") +
  xlab("# of genes") +
  ylab("# of cells") +
  ggtitle("Quality Control") +
  theme(plot.title = element_text(size = base_font_size, face = "bold", hjust = 0.5))

p2 <- plotRLE(
    x,
    exprs_mats = list(Raw = "logcounts_raw", 
                      CPM = "logcounts"),
    exprs_logged = c(TRUE, TRUE),
    colour_by = "total_features") +
  scale_color_continuous(low = "black", high = "black") +
  scale_fill_continuous(low = "black", high = "black") +
  guides(fill=FALSE, color=FALSE) +
  # guides(fill=guide_legend(title="# of genes")) + 
  # guides(color=guide_legend(title="# of genes")) + 
  theme_classic(base_size = base_font_size) +
  xlab("Cells") +
  ylab("RLE") +
  ggtitle("Normalisation") +
  theme(plot.title = element_text(size = base_font_size, face = "bold", hjust = 0.5))

p2$layers[[1]]$aes_params$colour <- "black"
p2$layers[[4]]$aes_params$size <- 0.1

p3 <- BrenneckeGetVariableGenes_ggplot(assays(x)[["norm"]], font_size = base_font_size) +
  ggtitle("Feature Selection") +
  theme(plot.title = element_text(size = base_font_size, face = "bold", hjust = 0.5))

x <- runPCA(x, ntop = 500)

p4 <- ggplot(as.data.frame(reducedDim(x)), 
       aes(x = PC1, 
           y = PC2)) + 
  geom_point() + 
  theme_classic(base_size=base_font_size) +
  guides(colour=FALSE) +
  ggtitle("Dimensionality Reduction") +
  theme(plot.title = element_text(size = base_font_size, face = "bold", hjust = 0.5))

dr <- x@reducedDims$PCA
d <- as.matrix(dist(dr))
p5 <- pheatmap(d, show_rownames = F, 
               show_colnames = F,
               treeheight_row = 20,
               treeheight_col = 20,
               fontsize = base_font_size-2,
               main = "Cell-Cell Distances")
p5 <- p5$gtable



set.seed(101)
c1 <- kmeans(dr, centers=5, nstart=100)

dat <- as.data.frame(reducedDim(x))
dat$Cluster <- factor(c1$cluster)
p6 <- ggplot(dat, 
       aes(x = PC1, 
           y = PC2, color = Cluster)) + geom_point() +
  scale_colour_brewer(palette = "Set3") +
  geom_point(data = as.data.frame(c1$centers), color = "black", size = 4, shape = 17) +
  theme_classic(base_size=base_font_size) + 
  guides(colour=FALSE) +
  ggtitle("Unsupervised Clustering") +
  theme(plot.title = element_text(size = base_font_size, face = "bold", hjust = 0.5))

arr <- ggdraw() + draw_image("arrow.png", scale = 0.9)

first_row <-  plot_grid(p1, arr, p2, arr, p3, arr, ncol = 6, rel_widths = c(1, 0.2, 1, 0.2, 1, 0.2))
second_row <- plot_grid(arr, p4, arr, p5, arr, p6, ncol = 6, rel_widths = c(0.2, 1, 0.2, 1, 0.2, 1))

plot_grid(NULL, first_row, NULL, second_row, NULL, ncol = 1, rel_heights = c(0.1, 1, 0.2, 1, 0.1))
ggsave("pdf/fig1.pdf", w = 9, h = 6)
ggsave("png/fig1.png", w = 9, h = 6)
