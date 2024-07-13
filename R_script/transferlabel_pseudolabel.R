unannoted = subset(new_merge,SampleID == 'P1019S2')
unannoted = NormalizeData(unannoted)
new_merge = subset(new_merge,SampleID %in% c('P1013S2','P1015S2','P1018S1','P1020S1','RX01'))
new_merge = NormalizeData(new_merge) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
  
unannoted.anchor <- 
  FindTransferAnchors(reference = new_merge, query = unannoted, dims = 1:30,reference.reduction = "pca")
predictions_major <- TransferData(anchorset = unannoted.anchor, refdata = new_merge$cl_major, dims = 1:30)
predictions_major[,c(2:17)]= NULL
colnames(predictions_major) = 'cl_major'
predictions_minor <- TransferData(anchorset = unannoted.anchor, refdata = new_merge$cl_minor, dims = 1:30)
predictions_minor[,c(2:31)] = NULL
colnames(predictions_minor) = 'cl_minor'

unannoted <- AddMetaData(unannoted, metadata = predictions_major)
unannoted <- AddMetaData(unannoted, metadata = predictions_minor)

table(new$cl_major,new$SampleID)
table(new$cl_minor,new$SampleID)

new_anno$pseudo_pt = if_else(new_anno$SampleID %in% c('P1013S2','P1018S1'),'P01',
                                  if_else(new_anno$SampleID %in% c('P1015S2','P1020S1'),'P02',
                                          'P03'))
new_anno$pseudo_tp = if_else(new_anno$SampleID %in% c('P1013S2','P1015S2','P1019S2'),'post','pre')

new_anno$pseudo_sample = paste0(new_anno$pseudo_pt,'_',new_anno$pseudo_tp)

sceasy::convertFormat(new_anno, from = 'seurat',to  ='anndata',  out = '/home/zhepan/new_anno.h5ad',main_layer = 'counts')
