library(SoupX)

path = '~/Project/MultiOmics/data/snRNA/Result/P1018S1'
dir.create(file.path(path,'soupx_output'))
setwd(path)

sc = load10X(path)

seu = CreateSeuratObject(sc$toc) %>% SCTransform(.) %>% RunPCA(., npcs = 40) %>%
  FindNeighbors(.,dims = 1:30) %>% FindClusters(., resolution = 0.8) %>% 
  RunUMAP(.,dims = 1:30)

matx = seu@meta.data

sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)
out = adjustCounts(sc,roundToInt = T)

new = CreateSeuratObject(out)


saveRDS(new,file = 'soupx_output/soupx.rds')

source("~/Project/MultiOmics/code/func/convertSeu5Format.R")
convertSeu5Format(new, savepath = 'soupx_output/soupx.h5ad')
