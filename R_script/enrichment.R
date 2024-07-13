#enrichment
seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

seu_sub = seu %>% subset(.,cl_major == 'Keratinocyte')
count = seu@assays$RNA@counts

library(GSVA)
library(GSEABase)

ncores = 4

gmt_file = '~/Reference/gmt/c5.go.v2023.2.Hs.symbols.gmt'
geneSets <- getGmt(gmt_file) 

gene_select = c(names(geneSets) %>% grep('COMPLEMENT',.,value = T), 
                names(geneSets) %>% grep('IMMUN',.,value = T),
                names(geneSets) %>% grep('METABOLISM',.,value = T))

geneSets = geneSets[gene_select]
gsva_es <- gsva(as.matrix(count), min.sz = 3,geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) 
signature_exp<-as.matrix(gsva_es) %>% t()
saveRDS(signature_exp, file = 'myeloid_meta_score.rds')

# seu = subset macrophage cluster
seu$cl_minor = gsub('Macrophage-','',seu$cl_minor)
seu$cl_minor = paste0(seu$cl_minor,'-Mac')

p1 = VlnPlot(seu, features = 'GOBP_COMPLEMENT_ACTIVATION',group.by = 'cl_minor',pt.size = 0)+
  labs(title = 'Complement Activation', caption = '', x = '')+
  ggsci::scale_fill_nejm()

p2 = VlnPlot(seu, features = 'GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION',group.by = 'cl_minor',pt.size = 0)+
  labs(title = 'Regulation of\n Complement Activation', caption = '', x = '')+
  ggsci::scale_fill_nejm()

p3 = 
  VlnPlot(seu, features = 'GOBP_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE',group.by = 'cl_minor',pt.size = 0)+
  labs(title = 'Suppressive\n Immune Response', caption = '', x = '')+
  ggsci::scale_fill_nejm()

p = (p1|p2|p3) + plot_layout(widths = c(1,1,1), guides = 'collect')
ggsave('ssgsea_mac.pdf',p,width = 11.38,height = 3.16,dpi = 300)
# Msigdb:C7 immune 

# metabolism
# scMetabolism

# run deg first
source("~/Project/MultiOmics/code/func/rungsea.R")


res = markers %>% run.gsea(.,gmt = '/home/zhepan/Reference/gmt/KEGG_metabolism_nc.gmt', 
                           clusters ='post', method = 'gsea' )

df = res$res@result

library(GseaVis)

sigs = c(
  'Oxidative phosphorylation','Fatty acid elongation','Tryptophan metabolism',
  'Arginine and proline metabolism','Propanoate metabolism','Glutathione metabolism')

for(i in sigs){
  p = gseaNb(object = res$res,geneSetID = i,addPval = F,subPlot = 2)
  ggsave(paste0(i,'_gsea.png'),p, width = 4.98, height = 4.26)
}

