#enrichment
seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

seu_sub = seu %>% subset(.,cl_major == 'Keratinocyte')
count = seu_sub@assays$RNA@counts

library(GSVA)
library(GSEABase)

ncores = 4

gmt_file = '~/Reference/gmt/KEGG_metabolism_nc.gmt'

geneSets <- getGmt(gmt_file) 
gsva_es <- gsva(as.matrix(count), min.sz = 3,geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) 
signature_exp<-as.matrix(gsva_es) %>% t()
saveRDS(signature_exp, file = 'meta_score.rds')

# Msigdb:C7 immune 


# metabolism
# scMetabolism