markers <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_cl_minor_marker.rds")
seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")
source("~/Project/MultiOmics/code/func/dotplot_config.R")
seu = seu %>% NormalizeData() %>% ScaleData()

Idents(seu) = seu$cl_minor
markers = FindAllMarkers(seu)

topmarker = markers %>% group_by(cluster) %>% top_n(20,pct.1-pct.2)

select = 'Monocyte-CXCL8'
topmarker %>% dplyr::filter(.,cluster == select)

saveRDS(markers,file = 'markers.rds')

# TREM2-Mac
marker_genes_dict = list(
  "Macrophage-MAML2" = c("MAML2", "LRMDA",'XIST'),
  "Macrophage-TREM2" = c("TREM2","C1QC","MS4A4A",'CSF1R','SPP1'),
  "Macrophage-IGFBP7" = c("IGFBP7"),
  "Monocyte-FCN1" = c("S100A8", "FCN1"),
  "Monocyte-CXCL8" = c("CXCL3","CXCL8")
)

marker_genes_dict = list(
  "Macrophage-TREM2" = c('SLC7A7','SLC7A6','SLC7A1','ACY1','GOT2','GOT1','ARG2','ASL')
)

dotplot_cl_order = c('Macrophage-TREM2','Macrophage-IGFBP7','Macrophage-MAML2')

dotplot_cl_order = c('Macrophage-TREM2','Macrophage-IGFBP7','Macrophage-MAML2','Monocyte-FCN1','Monocyte-CXCL8')

marker_genes_dict = list(
  "KC-Basal" = c("KRT15", "COL17A1"),
  "KC-Suprabasal" = c("KRT10","KRT1","CALN1"),
  "KC-Spinous" = c("EPS8L1","SPINK5","IL18"),
  "KC-Proliferating" = c("DIAPH3", "ASPM"),
  "KC-Cycling" = c("KRT6A","S100A8"),
  "Mast cell" = c("HDC","CPA3"),
  "Melanocyte" = c("TRPM1", "DCT"),
  "Endo-Vas" = c("SELE","ADGRL4"), 
  "Endo-Lymph" = c('VWF','PECAM1',"PKHD1L1", "CCL21"), 
  "PVL" = c("ACTA2"),
  "Fibroblast" = c("PDGFRA","COL1A1", "COL3A1"), 
  "Schwann" = c("NRXN1", "CADM2"),
  "Eccrine gland" = c("EDAR"),
  "DC" = c("CD1C", "CCL22", "CSF2RA"),  
  "Macrophage" = c("CD163","MS4A7","SPP1","FOLR2"),
  "LC" = c("PRKCB","CD207"),
  "NK" = c("GNLY","CD96","KLRC1"),
  "T-CD8" = c("GZMA","GZMK","IFNG","SCML4"),
  "T-CD4" = c('PTPRC','CD2','CD3G',"ITK","ICOS","CD40LG")
)

dotplot_cl_order = c('T-CD4','T-CD8','NK','Macrophage','DC','LC','Mast cell',
                     'KC-Basal','KC-Suprabasal','KC-Spinous','KC-Proliferating','KC-Cycling',
                     'Eccrine gland','Melanocyte',
                     'Endo','Endo-Lymph','Endo-Vas','PVL','Fibroblast',
                     'Schwann')

Idents(seu) = factor(seu$myeloid_minor_order, levels = dotplot_cl_order)
levels(marker_genes_dict) = dotplot_cl_order

dotplot_gene = rev(marker_genes_dict[levels(marker_genes_dict)]) %>% unlist() %>% unname()



p = DotPlot(seu, features = dotplot_gene)

pct_threshold = 0
exp_threshold = 0

p1 = 
  ggplot(data = p$data %>% dplyr::filter(.,pct.exp>pct_threshold & avg.exp.scaled>exp_threshold),
       aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
  geom_point()+
  coord_flip()+
  scale_color_viridis_c()+
  p$theme+
  theme(plot.background = element_rect(fill = "white"),
        plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
        panel.border = element_rect(fill = NA, color = "black", size = 0.7),
        aspect.ratio = 2.4,
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size = 8,angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 8)) + 
  labs(x = '', y = '')+
  guides(
  fill = guide_legend(title=""), 
  size = guide_legend(title='Percentage\nExpressed',override.aes = list(size=5)),
  colour = guide_colorbar(title='Scaled\nexpression', override.aes = list(size=5),frame.colour = "black",label = TRUE,label.vjust = 1,ticks = FALSE)
)
p1
ggsave('mac_dotplot.pdf',p1, width = 2.74,height = 4.62,scale = 1,dpi = 300)

# subset
seu <- readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

sub_cl = 'T cell'
sub_seu = seu %>% subset(.,cl_major == sub_cl) %>% NormalizeData() %>% ScaleData()

top_m_list[[sub_cl]] %>% dplyr::filter(.,cluster == 'NK') %>% print(n = 20)

p1 = plotDot(sub_seu, marker_gene_dict = t_marker_genes_dict, dotplot_cl_order = t_dotplot_cl_order)

ggsave(paste0(sub_cl,'_dotplot.png'),p1, width = 8.15,height = 3.9,scale = 1,dpi = 300)


plotDot = function(seu, marker_gene_dict, dotplot_cl_order, pct_threshold = 5 , exp_threshold = 0.5, aspect_ratio = 0.5){
  
  Idents(seu) = factor(seu$cl_minor, levels = rev(dotplot_cl_order))
  levels(marker_gene_dict) = dotplot_cl_order
  
  dotplot_gene = marker_gene_dict[levels(marker_gene_dict)] %>% unlist() %>% unname()
  
  p = DotPlot(seu, features = dotplot_gene)
  
  p1 = 
    ggplot(data = p$data %>% dplyr::filter(.,pct.exp>pct_threshold & avg.exp.scaled>exp_threshold),
           aes(x = features.plot, y = id, color = avg.exp.scaled, size = pct.exp))+
    geom_point()+
    coord_flip()+
    scale_color_viridis_c()+
    p$theme+
    theme(plot.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.25,0,0.25,1), "cm"), 
          panel.border = element_rect(fill = NA, color = "black", size = 0.7),
          aspect.ratio = aspect_ratio,
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          axis.ticks = element_blank(),
          axis.text.x = element_text(size = 8,angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 8),
          legend.position = 'right') + 
    labs(x = '', y = '')+
    guides(
      fill = guide_legend(title=""), 
      size = guide_legend(title='Percentage\nExpressed'),
      colour = guide_colorbar(title='Scaled\nexpression', override.aes = list(size=5),frame.colour = "black",label = TRUE,label.vjust = 1,ticks = FALSE)
    )
  return(p1)
  
}
plotDot(seu, marker_genes_dict, dotplot_cl_order, pct_threshold = 0 , exp_threshold = 0.5, aspect_ratio = 2)
  