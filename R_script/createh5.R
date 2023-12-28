library(scCustomize)

paths = paste0('~/Project/MultiOmics/data/skin/',c('SKIN-A1002','SKIN-A1007','SKIN-A1009',
                                                   'SKIN-B1002','SKIN-B1007','SKIN-B1009'))
for (path in paths) {
  
  setwd(path)

  Create_10X_H5(raw_data_file_path = "filtered_feature_bc_matrix", save_file_path = getwd(), save_name = "filtered_feature_bc_matrix")
  Create_10X_H5(raw_data_file_path = "raw_feature_bc_matrix", save_file_path = getwd(), save_name = "raw_feature_bc_matrix")
}