source('/home/zhepan/Project/MultiOmics/code/func/obj_create.R')

#for(i in c('cellbender','cellranger','soupx')){
#obj_create(
#  type = i,
#  datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result',
#  savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object',
#  PatientID = 'P1013',
#  NeoChemoRes = 'PCR',
#  NeoRadRes = 'PCR',
#  SampleID = 'P1013S2',
#  SampleTimepoint = 'S2',
#  SampleMethod = 'surgery', # biopsy
#  SampleDate = '20231102',
#  Kit = 'sc5'
#)
#}

for (i in c( 'cellbender')) {
  obj_create(
    project = 'skin',
    type = i,
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1002',
    SampleID = 'SKIN-A1002',
    SampleType = 'N',
    SampleDate = '20230111',
    Kit = 'sn3'
  )

}
  obj_create(
    project = 'skin',
    type = 'cellbender',
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1002',
    SampleID = 'SKIN-B1002',
    SampleType = 'H',
    SampleDate = '20230111',
    Kit = 'sn3'
  )
  
  obj_create(
    project = 'skin',
    type = 'cellbender',
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1007',
    SampleID = 'SKIN-A1007',
    SampleType = 'N',
    SampleDate = '20230808',
    Kit = 'sn3'
  )
  
  obj_create(
    project = 'skin',
    type = 'cellbender',
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1007',
    SampleID = 'SKIN-B1007',
    SampleType = 'H',
    SampleDate = '20230808',
    Kit = 'sn3'
  )
  
  
  obj_create(
    project = 'skin',
    type = 'cellbender',
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1009',
    SampleID = 'SKIN-A1009',
    SampleType = 'N',
    SampleDate = '20230726',
    Kit = 'sn3'
  )
  
  obj_create(
    project = 'skin',
    type = 'cellbender',
    datapath = '/home/zhepan/Project/MultiOmics/data/skin',
    savepath = '/home/zhepan/Project/MultiOmics/data/skin',
    PatientID = 'P1009',
    SampleID = 'SKIN-B1009',
    SampleType = 'H',
    SampleDate = '20230726',
    Kit = 'sn3'
  )


# test for decontamination effect -----------------

seu = subset(seu, nCount_RNA <20000 & nCount_RNA >500 & nFeature_RNA<6000 & scDblFinder.class == 'singlet')
seu = seu %>% NormalizeData(.) %>% FindVariableFeatures(.) %>% ScaleData(.) %>% RunPCA() %>%  RunUMAP(.,dims = 1:30)

#library(viridis)
p4 = FeaturePlot(seu,c('KRT14','KRT10','KRT1'), ncol = 3,pt.size = 0.1)&scale_color_viridis_c()
p0/p1/p2
