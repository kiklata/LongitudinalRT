library(SoupX)

paths = paste0('~/Project/MultiOmics/data/skin/',c('SKIN-A1002','SKIN-A1007','SKIN-A1009',
                                                   'SKIN-B1002','SKIN-B1007','SKIN-B1009'))
for (path in paths) {
  
setwd(path)

sc = load10X(paths[1])

seu = CreateSeuratObject(sc$toc) %>% SCTransform(.) %>% RunPCA(., npcs = 40) %>%
  FindNeighbors(.,dims = 1:30) %>% FindClusters(., resolution = 0.8) %>% 
  RunUMAP(.,dims = 1:30)

matx = seu@meta.data

sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
sc = autoEstCont(sc)
out = adjustCounts(sc,roundToInt = T)

new = CreateSeuratObject(out)


saveRDS(new,file = 'raw/soupx.rds')

source("~/Project/MultiOmics/code/func/convertSeu5Format.R")
convertSeu5Format(new, savepath = 'raw/soupx.h5ad')
}