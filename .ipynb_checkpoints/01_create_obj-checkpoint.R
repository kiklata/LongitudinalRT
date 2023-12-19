source('/home/zhepan/Project/MultiOmics/code/func/obj_create.R')

obj_create(
  datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result',
  savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object',
  PatientID = 'P1013',
  NeoChemoRes = 'pCR',
  NeoRadRes = 'pCR',
  SampleID = 'P1013S2',
  SampleTimepoint = 'S2',
  SampleMethod = 'surgery',
  SampleDate = '20231101',
  Kit = 'sc5'
)
