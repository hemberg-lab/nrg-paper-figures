library(SingleCellExperiment)
library(scater)
library(cowplot)
library(statmod)
library(ggthemes)
library(igraph)
library(cccd)
library(Matrix)
library(GGally)
library(network)
library(sna)
library(intergraph)
library(RColorBrewer)
library(dendextend)
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

x <- runPCA(x, ntop = 500)

p1 <- ggplot(as.data.frame(reducedDim(x)), 
             aes(x = PC1, 
                 y = PC2,
                 color = x$cell_type2)) + 
  geom_point() + 
  theme_classic(base_size=base_font_size) +
  scale_color_tableau("colorblind10") +
  guides(colour=FALSE)

dr <- x@reducedDims$PCA
set.seed(101)
c1 <- kmeans(dr, centers=5, nstart=100)

dat <- as.data.frame(reducedDim(x))
dat$Cluster <- factor(c1$cluster)
p2 <- ggplot(dat, 
             aes(x = PC1, 
                 y = PC2, color = Cluster)) + geom_point() +
  scale_colour_brewer(palette = "Set3") +
  geom_point(data = as.data.frame(c1$centers), color = "black", size = 4, shape = 17) +
  theme_classic(base_size=base_font_size) + 
  guides(colour=FALSE)

set.seed(281)
kNN <- nng(x = dr, k=5)
adj_knn = get.adjacency(kNN)
snn <- adj_knn%*%t(adj_knn)
diag(snn) <- 0
sNN <- graph.adjacency(snn, mode="undirected")
sNN <- simplify(sNN)
louv <- igraph::cluster_louvain(sNN)
p31 <- ggnet2(sNN, mode = dr, node.size = 1.5, node.color = brewer.pal(9,"Set3")[louv$membership], edge.size = 0.25, edge.color = "black")

set.seed(281)
kNN <- nng(x = dr, k=10)
adj_knn = get.adjacency(kNN)
snn <- adj_knn%*%t(adj_knn)
diag(snn) <- 0
sNN <- graph.adjacency(snn, mode="undirected")
sNN <- simplify(sNN)
louv <- igraph::cluster_louvain(sNN)
p32 <- ggnet2(sNN, mode = dr, node.size = 1.5, node.color = brewer.pal(9,"Set3")[louv$membership], edge.size = 0.25, edge.color = "black")


hc <- dend <- dr %>% scale %>% dist %>% hclust
ord <- order(hc)
tr <- cutree(hc, k = 5)

g <- ggplot_build(p1)
dend <- dr %>% scale %>% dist %>% 
  hclust %>% as.dendrogram %>%
  set("leaves_pch", 19) %>%
  set("leaves_cex", 2) %>%
  set("leaves_col", brewer.pal(6,"Set3")[tr][ord]) %>% set("branches_lwd", 0.5) %>%
  set("labels", rep("", 268))

ggd1 <- as.ggdend(dend)
p4 <- ggplot(ggd1) 

first_row <-  plot_grid(p1, p2, p31, ncol = 3, labels = c("a", "b", "c"), rel_widths = c(1, 1, 1))
second_row <- plot_grid(p4, p32, ncol = 2, labels = c("d", "e"), rel_widths = c(2, 1))
plot_grid(first_row, second_row, ncol = 1)
ggsave("pdf/fig3.pdf", w = 9, h = 6)
ggsave("png/fig3.png", w = 9, h = 6)
