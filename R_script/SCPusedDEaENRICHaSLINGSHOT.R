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
