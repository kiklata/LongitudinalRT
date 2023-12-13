source('/home/zhepan/Project/MultiOmics/code/func/obj_create.R')

for(i in c('cellbender','cellranger','soupx')){
obj_create(
  type = i,
  datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result',
  savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object',
  PatientID = 'P1013',
  NeoChemoRes = 'PCR',
  NeoRadRes = 'PCR',
  SampleID = 'P1013S2',
  SampleTimepoint = 'S2',
  SampleMethod = 'surgery', # biopsy
  SampleDate = '20231102',
  Kit = 'sc5'
)
}

for(i in c('cellbender','cellranger','soupx')){
  obj_create(
    type = i,
    datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result',
    savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object',
    PatientID = 'P1015',
    NeoChemoRes = 'PCR',
    NeoRadRes = 'PCR',
    SampleID = 'P1015S2',
    SampleTimepoint = 'S2',
    SampleMethod = 'surgery', # biopsy
    SampleDate = '20231114',
    Kit = 'sc5'
  )
}

for(i in c('cellbender','cellranger','soupx')){
  obj_create(
    type = i,
    datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result',
    savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object',
    PatientID = 'P1018',
    NeoChemoRes = 'PCR',
    NeoRadRes = 'PCR',
    SampleID = 'P1018S1',
    SampleTimepoint = 'S1',
    SampleMethod = 'biopsy', # biopsy
    SampleDate = '20231107',
    Kit = 'sc5'
  )
}