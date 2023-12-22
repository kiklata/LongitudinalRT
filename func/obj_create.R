obj_create = function(type,
                      datapath,
                      savepath,
                      PatientID = NULL,
                      NeoChemoRes = NULL,
                      NeoRadRes = NULL,
                      SampleID = NULL,
                      SampleType = NULL,
                      SampleTimepoint = NULL,
                      SampleMethod = NULL,
                      SampleDate = NULL,
                      Kit = NULL) {
  

# only for test -----------------------------------------------------------
  #datapath = '/home/zhepan/Project/MultiOmics/data/snRNA/Result'
  #savepath = '/home/zhepan/Project/MultiOmics/data/snRNA/Object'
  #PatientID = 'P1015'
  #NeoChemoRes = 'NOTK'
  #NeoRadRes = 'NOTK'
  #SampleID = 'P1015S2'
  #SampleTimepoint = 'S2'
  #SampleMethod = 'surgery' # biopsy
  #SampleDate = '20231114'
  #Kit = 'sc5'
  
  
  # load pkgs ---------------------------------------------------------------
  
  suppressMessages(library(Seurat))
  suppressMessages(library(dplyr))
  suppressMessages(library(scDblFinder))
  
  
  # load data ---------------------------------------------------------------
  if(type =='cellranger'){
    count = Read10X_h5(file.path(datapath, SampleID, 'filtered_feature_bc_matrix.h5'))  
    seu = CreateSeuratObject(
      count,
      project = SampleID,
      min.cells = 3,
      min.features = 200
      )
  }else if(type == 'soupx'){
    seu = readRDS(file.path(datapath, SampleID, 'soupx_output','soupx.rds'))
  }else if(type == 'cellbender'){
    count = scCustomize::Read_CellBender_h5_Mat(file.path(datapath, SampleID, 'cellbender_output','cellbender_feature_bc_matrix_filtered.h5'))
    seu = CreateSeuratObject(
      count,
      project = SampleID,
      min.cells = 3,
      min.features = 200
    )
  }else if(type == 'skin'){
    count = Read10X(file.path(datapath, SampleID))
    seu = CreateSeuratObject(
      count,
      project = SampleID,
      min.cells = 3,
      min.features = 200
    )
  }
  
  if(type != 'skin'){
  seu$PatientID = PatientID
  seu$NeoChemoRes = NeoChemoRes
  seu$NeoRadRes = NeoRadRes
  seu$SampleID = SampleID
  seu$SampleTimepoint = SampleTimepoint
  seu$SampleMethod = SampleMethod
  seu$SampleDate = SampleDate
  seu$Kit = Kit
  }else {
    seu$PatientID = PatientID
    seu$SampleID = SampleID
    seu$SampleType = SampleType
    seu$SampleDate = SampleDate
    seu$Kit = Kit
  }
  
  # doublet detect ----------------------------------------------------------
  
  set.seed(42)
  sce = scDblFinder(as.SingleCellExperiment(seu))
  seu$scDblFinder.score = sce$scDblFinder.score
  seu$scDblFinder.class = sce$scDblFinder.class
  
  
  # some qc info ------------------------------------------------------------
  
  seu$percent_mt = PercentageFeatureSet(seu, pattern = '^MT-')
  seu$percent_hb=PercentageFeatureSet(seu, "^HB[^(P)]")
  
  cyclegenes = read.delim('/home/zhepan/Reference/regev_lab_cell_cycle_genes.txt')
  seu <- NormalizeData(seu)
  
  seu = CellCycleScoring(seu, s.features = cyclegenes[1:42, 1], g2m.features = cyclegenes[43:96, 1])
  
  #VlnPlot(seu,features = c('nCount_RNA','nFeature_RNA','percent_mt','S.Score'), group.by = 'scDblFinder.class', pt.size = 0)
  
  
  # convert to h5ad ---------------------------------------------------------
  dir.create(path = file.path(savepath, SampleID,'raw'), recursive = T)
  source("~/Project/MultiOmics/code/func/convertSeu5Format.R")  

  if(type == 'cellranger'){
    h5names = 'cellranger_doublet.h5ad'
    rdsnames = 'cellranger_doublet.rds'
  }else if(type == 'soupx'){
    h5names = 'soupx_doublet.h5ad'
    rdsnames = 'soupx_doublet.rds'
  }else if(type == 'cellbender'){
    h5names = 'cellbender_doublet.h5ad'
    rdsnames = 'cellbender_doublet.rds'
  }else if(type == 'skin'){
    h5names = 'cellranger_doublet.h5ad'
    rdsnames = 'cellranger_doublet.rds'
  }
  seu = DietSeurat(seu, layers = 'counts')
  convertSeu5Format(seu, savepaths = file.path(savepath, SampleID, 'raw',h5names))
  saveRDS(seu, file = file.path(savepath, SampleID, 'raw',rdsnames))

}

