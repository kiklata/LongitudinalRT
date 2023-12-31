library(celda)
library(scDblFinder)

paths = paste0('~/Project/MultiOmics/data/skin/',c('SKIN-A1002','SKIN-A1007','SKIN-A1009',
                                                   'SKIN-B1002','SKIN-B1007','SKIN-B1009'))
for (path in paths) {
  
  setwd(path)
  filter_count = Read10X('filtered_feature_bc_matrix')
  raw_count = Read10X('raw_feature_bc_matrix')
  
  decon = decontX(filter_count, background = raw_count)
  decon_count = decon$decontXcounts %>% round(.,0)
  seu = CreateSeuratObject(decon_count,min.cells = 3,min.features = 200)
  saveRDS(seu,file = 'decontx.rds')
}  
