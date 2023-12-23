library(celda)
library(scDblFinder)

count = cellranger_doublet@assays$RNA$counts

decon <- decontX(count)

new_count = decon$decontXcounts

seu = CreateSeuratObject(new_count %>% round(.,0))
sce = scDblFinder(as.SingleCellExperiment(seu))
seu$scDblFinder.score = sce$scDblFinder.score
seu$scDblFinder.class = sce$scDblFinder.class

seu = subset(seu, nCount_RNA <20000 & nCount_RNA >500 & nFeature_RNA<6000 & scDblFinder.class == 'singlet')
seu = seu %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA() %>%  RunUMAP(.,dims = 1:30)

old = subset(cellranger_doublet, nCount_RNA <20000 & nCount_RNA >500 & nFeature_RNA<6000 & scDblFinder.class == 'singlet')
old = old %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA() %>%  RunUMAP(.,dims = 1:30)

library(viridis)
p = (FeaturePlot(seu,c('KRT14','KRT10','KRT1'), ncol = 3,pt.size = 0.1)/FeaturePlot(old,c('KRT14','KRT10','KRT1'), ncol = 3,pt.size = 0.1))&scale_color_viridis_c()
ggsave('/home/zhepan/Project/MultiOmics/code/figures/decontx_without_raw_A1002.png',p,width = 10,height = 8)
