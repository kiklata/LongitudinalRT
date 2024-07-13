obj_create = function(project,
                      type,
                      datapath,
                      savepath,
                      PatientID = NULL,
                      SampleID = NULL,
                      SampleType = NULL,
                      SampleTimepoint = NULL,
                      SampleMethod = NULL,
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
  
  
  if (project == 'tumor') {
    cellranger_path = file.path(datapath, SampleID, 'filtered_feature_bc_matrix.h5')
    soupx_path = file.path(datapath, SampleID, 'soupx.rds')
    decontx_path = file.path(datapath, SampleID, 'decontx.rds')
    cellbender_path = file.path(
      datapath,
      SampleID,
      'cellbender_output',
      'cellbender_feature_bc_matrix_filtered.h5'
    )
    
    if (type == 'cellranger') {
      count = Read10X_h5(cellranger_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'soupx') {
      seu = readRDS(soupx_path)
    } else if (type == 'cellbender') {
      count = scCustomize::Read_CellBender_h5_Mat(cellbender_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'decontx') {
      seu = readRDS(decontx_path)
    }
    seu$PatientID = PatientID
    seu$SampleID = SampleID
    seu$SampleTimepoint = SampleTimepoint
    seu$SampleMethod = SampleMethod
    seu$Kit = Kit
    
    
  } else if (project == 'skin') {
    cellranger_path = file.path(datapath, SampleID,'raw/filtered_feature_bc_matrix')
    soupx_path = file.path(datapath, SampleID, 'raw/soupx.rds')
    cellbender_path = file.path(datapath,
                                SampleID,
                                'raw/cellbender_filtered.h5')
    decontx_path = file.path(datapath, SampleID, 'raw/decontx.rds')
    
    if (type == 'cellranger') {
      count = Read10X(cellranger_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'soupx') {
      seu = readRDS(soupx_path)
    } else if (type == 'cellbender') {
      count = scCustomize::Read_CellBender_h5_Mat(cellbender_path)
      seu = CreateSeuratObject(
        count,
        project = SampleID,
        min.cells = 3,
        min.features = 200
      )
    } else if (type == 'decontx') {
      seu = readRDS(decontx_path)
    }
    seu$PatientID = PatientID
    seu$SampleID = SampleID
    seu$SampleType = SampleType
    seu$SampleDate = SampleDate
    seu$Kit = Kit
  } else{
    stop("please specify project")
  }
  
  # doublet detect ----------------------------------------------------------
  
  set.seed(42)
  sce = scDblFinder(as.SingleCellExperiment(seu))
  seu$scDblFinder.score = sce$scDblFinder.score
  seu$scDblFinder.class = sce$scDblFinder.class
  
  
  # some qc info ------------------------------------------------------------
  
  seu$percent_mt = PercentageFeatureSet(seu, pattern = '^MT-')
  seu$percent_hb = PercentageFeatureSet(seu, pattern = "^HB[^(P)]")
  seu$percent_rb = PercentageFeatureSet(seu, pattern = "^RP[SL]")
  
  cyclegenes = read.delim('/home/zhepan/Reference/regev_lab_cell_cycle_genes.txt',header = F)
  seu <- NormalizeData(seu)
  
  seu = CellCycleScoring(seu, s.features = cyclegenes[1:42, 1], g2m.features = cyclegenes[43:96, 1])
  
  #VlnPlot(seu,features = c('nCount_RNA','nFeature_RNA','percent_mt','S.Score'), group.by = 'scDblFinder.class', pt.size = 0)
  
  
  # convert to h5ad ---------------------------------------------------------
  if(!dir.exists(paths = file.path(savepath, SampleID,'raw'))){
    dir.create(path = file.path(savepath, SampleID,'raw'),recursive = T)}
  else {NULL}
  
  if (type == 'cellranger') {
    h5names = 'cellranger_doublet.h5ad'
    rdsnames = 'cellranger_doublet.rds'
  } else if (type == 'soupx') {
    h5names = 'soupx_doublet.h5ad'
    rdsnames = 'soupx_doublet.rds'
  } else if (type == 'cellbender') {
    h5names = 'cellbender_doublet.h5ad'
    rdsnames = 'cellbender_doublet.rds'
  } else if (type == 'decontx') {
    h5names = 'decontx_doublet.h5ad'
    rdsnames = 'decontx_doublet.rds'
  } 
  
  DietSeurat(seu,counts = TRUE,data = TRUE,scale.data = FALSE)
  sceasy::convertFormat(seu, from="seurat", to="anndata",outFile = file.path(savepath, SampleID, 'raw', h5names), drop_single_values = F, main_layer = 'counts')
  saveRDS(seu, file = file.path(savepath, SampleID, 'raw', rdsnames))
  
}
