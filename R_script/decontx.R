library(celda)

paths = paste0('/home/zhepan/Project/MultiOmics/data/tumor/scRNA/Result/',c('P1024S2'))

for (path in paths) {
  
  setwd(path)
  filter_count = Read10X('filtered_feature_bc_matrix')
  raw_count = Read10X('raw_feature_bc_matrix')
  
  decon = decontX(filter_count, background = raw_count)
  decon_count = decon$decontXcounts %>% round(.,0)
  seu = CreateSeuratObject(decon_count,min.cells = 3,min.features = 200)
  saveRDS(seu,file = 'decontx.rds')
}  
