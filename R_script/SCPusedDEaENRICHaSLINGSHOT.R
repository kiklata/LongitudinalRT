tcell_count <- readRDS("~/Project/MultiOmics/data/snRNA/Object/summary/Immune/tcell_count.rds")
tcell_cluster <- read.delim("~/Project/MultiOmics/data/snRNA/Object/summary/annotation/tcell_cluster.tsv", row.names=1)
tcell_count = AddMetaData(tcell_count,tcell_cluster)

all.count =tcell_count@assays$RNA@counts

mt.name = grep ("^MT-", rownames(all.count),value = T)
rp.name = grep("^RP[L|S]",rownames(all.count),value = T)
hb.name = grep("^HB[^(P)]",rownames(all.count),value = T)

new.count = all.count[!rownames(all.count) %in% hb.name,]

new_tcell = CreateSeuratObject(new.count, meta.data = tcell_count@meta.data)

new_tcell = NormalizeData(new_tcell,assay = 'RNA')

# DE -----------------------------------------------------------------
library(SCP)

cd4t_sub <- RunDEtest(srt = new_tcell %>% subset(.,Tcell_minor == 'CD4-Treg-FOXP3'), group_by = "SampleTimepoint",group1 = 'S2', 
                      fc.threshold = 1, only.pos = FALSE)
cd8t_sub <- RunDEtest(srt = new_tcell %>% subset(.,Tcell_minor == 'CD8-Trm-ZNF683'), group_by = "SampleTimepoint", group1 = 'S2', 
                      fc.threshold = 1, only.pos = FALSE)
voltheme = theme(
  plot.title = element_text(hjust = 0.5),
  strip.background = element_blank(),
  strip.text.x = element_blank()
  )

VolcanoPlot(srt = cd4t_sub)+labs(title = 'CD4-Treg-FOXP3')+voltheme
  
VolcanoPlot(srt = cd8t_sub)+labs(title = 'CD8-Trm-ZNF683')+voltheme


# enrichment --------------------------------------------------------------
cd8t_major <- RunDEtest(srt = new_tcell %>% subset(.,Tcell_major == 'CD8-T'), group_by = "SampleTimepoint", 
                      fc.threshold = 1, only.pos = FALSE)

deg = cd8t_major@tools$DEtest_SampleTimepoint$AllMarkers_wilcox

source("~/Project/MultiOmics/code/func/rungsea.R")

res = deg %>% run.gsea(.,cluster ='S2', method = 'gsea' )
df = res$res@result
df$group = if_else(df$NES>0, 'post','pre')

plotdf = df %>% dplyr::filter(.,abs(NES)>0.5)

ggplot(plotdf)+
  geom_col(aes(reorder(Description,NES),y = NES,fill=group))+
  scale_fill_manual(values = c("#8D4873","#1084A4"))+
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(color="black",size=10),
        axis.line.x = element_line(color='black'),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))+
  coord_flip()+
  geom_segment(aes(y=0, yend=0,x=0,xend=nrow(plotdf)+1))+
  geom_text(data = plotdf[plotdf$NES> 0.5,],aes(x=Description, y=-0.01, label=Description),
            hjust=1, size=2)+
  geom_text(data = plotdf[plotdf$NES< -0.5,],aes(x=Description, y=0.01, label=Description),
            hjust=0, size=2)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, .5, 1))+
  labs(title = 'CD8-T',x='', y='Normalized Enrichment Score')


# slingshot ---------------------------------------------------------------

new_tcell@reductions = tcell_cluster@reductions
cd4 = subset(new_tcell, Tcell_major == 'CD4-T')

tcd4cell_maj <- RunSlingshot(srt = cd4, group.by = "Tcell_minor", reduction = "UMAP")

FeatureDimPlot(tcd4cell_maj, features = paste0("Lineage", 1), reduction = "UMAP", theme_use = "theme_blank")
CellDimPlot(tcd4cell_maj, group.by = "Tcell_minor", reduction = "UMAP", lineages = paste0("Lineage", 1))+ labs(title = '',caption = '')
  +theme(legend.position = 'bottom')


tcd4cell_maj[["RNA3"]] <- as(object = tcd4cell_maj[["RNA"]], Class = "Assay")
DefaultAssay(tcd4cell_maj) <- "RNA3"
tcd4cell_maj[["RNA"]] <- NULL
tcd4cell_maj <- RenameAssays(object = tcd4cell_maj, RNA3 = 'RNA')

tcd4cell_maj <- RunDynamicFeatures(srt = tcd4cell_maj, lineages = c("Lineage1"), n_candidates = 200)

ht <- DynamicHeatmap(
  srt = tcd4cell_maj, lineages = c("Lineage1"),feature_split_by = 'Lineage1',
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1", 
  heatmap_palette = "viridis", cell_annotation = "Tcell_minor",
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot)

DynamicPlot(
  srt = tcd4cell_maj, lineages = c("Lineage1"), group.by = "Tcell_minor",
  features = c("FOXP3",'IL7R'),
  compare_lineages = TRUE, compare_features = FALSE
)


# myeloid -----------------------------------------------------------------

## scp
ht1 <- GroupHeatmap(pancreas_sub,
                    features = c(
                      "Sox9", "Anxa2", "Bicc1", # Ductal
                      "Neurog3", "Hes6", # EPs
                      "Fev", "Neurod1", # Pre-endocrine
                      "Rbp4", "Pyy", # Endocrine
                      "Ins1", "Gcg", "Sst", "Ghrl" # Beta, Alpha, Delta, Epsilon
                    ),
                    group.by = c("CellType", "SubCellType")
)

count =myeloid_count@assays$RNA@counts

mt.name = grep ("^MT-", rownames(count),value = T)
rp.name = grep("^RP[L|S]",rownames(count),value = T)
hb.name = grep("^HB[^(P)]",rownames(count),value = T)

new.count = count[!rownames(count) %in% hb.name,]

new_myeloid = CreateSeuratObject(new.count, meta.data = myeloid_count@meta.data)

new_myeloid = NormalizeData(new_myeloid,assay = 'RNA')
new_myeloid = ScaleData(new_myeloid)

Idents(new_myeloid) = new_myeloid$Myeloid_minor
all_marker = FindAllMarkers(new_myeloid)

tes = all_marker %>% dplyr::group_by(cluster) %>% dplyr::filter(p_val_adj <0.05) %>% 
  top_n(n = 5, wt = avg_log2FC)

removeg = HGNChelper::checkGeneSymbols(x = tes$gene) %>% dplyr::filter(Approved == T) 

tes_l = tes[tes$gene %in% removeg$x,]

vag_exp = AverageExpression(obj, 
                            assays = "RNA", 
                            features = tes_l$gene,
                            group.by = "Myeloid_minor",
                            layer = "data")

celltype = unique(FetchData(obj, vars = c("myeloid_major","Myeloid_minor")))

celltype$myeloid_major = factor(celltype$myeloid_major, 
                                levels = c("Macrophage",'Monocyte','DC'))

celltype = celltype[order(celltype$myeloid_major), ]

celltype$Myeloid_minor = factor(celltype$Myeloid_minor, 
                                levels = unique(celltype$Myeloid_minor))

tes_l$cluster = factor(tes_l$cluster, levels = levels(celltype$Myeloid_minor))
tes_l = tes_l[order(tes_l$cluster),]
tes_l$gene = factor(tes_l$gene, levels = tes_l$gene)

dat = base::apply(vag_exp$RNA, 1, function(x) (x - mean(x)) / sd(x)) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene') %>% 
  reshape2::melt()

dat$Gene = factor(dat$Gene, levels = rev(levels(tes_l$gene)))

dat$variable = factor(dat$variable, levels = levels(celltype$Myeloid_minor))

dat = dat[order(dat$Gene), ]

heatmap_color = RColorBrewer::brewer.pal(name = "RdBu",n = 11)

pal = rev(colorRampPalette(heatmap_color)(500))
label1 = levels(celltype$Myeloid_minor)[seq(1, length(levels(celltype$Myeloid_minor)), by = 2)]
label2 = levels(celltype$Myeloid_minor)[seq(2, length(levels(celltype$Myeloid_minor)), by = 2)]


p1 = ggplot(dat, aes(as.numeric(variable), Gene, fill=value))+
  geom_tile() +
  scale_y_discrete(expand = c(0, 0))+
  scale_x_continuous(expand = c(0, 0),
                     breaks = seq(1, length(levels(celltype$Myeloid_minor)), by = 2),
                     labels = label1,
                     sec.axis = dup_axis(
                       breaks = seq(2, length(levels(celltype$Myeloid_minor)), by = 2),
                       labels = label2)
  ) +
scale_fill_gradientn(colors = pal, limits = c(-2.5, 3), name = "Z Score") +
  geom_hline(yintercept = as.numeric(cumsum(rev(table(tes_l$cluster)[-1])) + .5), linetype = 2)+
  geom_vline(xintercept = as.numeric(cumsum(table(celltype$myeloid_major)) + .5), linetype = 2)+
  theme(text = element_text(face = "bold"),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y=element_blank(),
        axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = .5),
        axis.text.x.bottom = element_text(angle = 90, hjust = 1, vjust = .5)
  )

tes_l$x = rep(c(1,2,3),31)[1:92]
tes_l$y = tes_l$x
initindex = 0
for ( i in names(table(tes_l$cluster))) {
  for (k in 1:nrow(tes_l[tes_l$cluster == i,])) {
    tes_l[tes_l$cluster == i,][k,'y'] = 31-initindex
    if( tes_l[tes_l$cluster == i,][k,'x'] == 3){
      initindex = initindex+1
    }else{NULL}
  }
}


p2 = ggplot(tes_l, aes(x,y,fill = cluster))+
  geom_tile()+
  geom_text(aes(label = gene), family = "serif", fontface = "italic") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_brewer(palette = "Pastel2") +
  theme(text = element_text(),
        axis.text = element_blank(),
        axis.title =element_blank(),
        axis.ticks = element_blank(),
        legend.position="none")+
  scale_x_continuous(expand = c(0,0)) + 
  geom_hline(yintercept = as.numeric(cumsum(rev(table(tes_l$cluster)[-1]/3))) + .5, color = "white")

## markergene number need manula select

library(patchwork)

pic.heatmap = p2 + p1 + plot_layout(ncol = 2, widths  = c(1, 3))

pdf("test.pdf",w=10.5,h=10)
pic.heatmap & theme(plot.margin = margin(0,0,0,0))
dev.off()