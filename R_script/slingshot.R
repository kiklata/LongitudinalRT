library(slingshot)

sds <- slingshot(Embeddings(seu, "umap"), clusterLabels = seu$seurat_clusters, 
                 start.clus = 4, stretch = 0)
library(Seurat) 
library(SingleCellExperiment)
sce <- as.SingleCellExperiment(seuratObject)
sce <- slingshot(sce, clusterLabels = ident, reducedDim = "PCA",
                 allow.breaks = FALSE)

library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

lin1 <- getLineages(rd, cl, start.clus = '1')
crv1 <- getCurves(lin1)
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')


library(Polychrome)
library(ggbeeswarm)
library(ggthemes)

# this define the cluster color. You can change it with different color scheme.
my_color <- createPalette(length(levels(sce$ident)), c("#010101", "#ff0000"), M=1000)
names(my_color) <- unique(as.character(sce$ident))

slingshot_df <- data.frame(colData(sce))

# re-order y-axis for better figure: This should be tailored with your own cluster names
# slingshot_df$ident = factor(slingshot_df$ident, levels=c(4,2,1,0,3,5,6))

ggplot(slingshot_df, aes(x = slingPseudotime_1, y = ident, 
                         colour = ident)) +
  geom_quasirandom(groupOnX = FALSE) + theme_classic() +
  xlab("First Slingshot pseudotime") + ylab("cell type") +
  ggtitle("Cells ordered by Slingshot pseudotime")+scale_colour_manual(values = my_color)
