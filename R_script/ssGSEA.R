## ssGSEA

curated_meta = read.delim('~/Project/MultiOmics/data/snRNA/Object/summary/annotation/cellranger_filter_anno.tsv',row.names = 1)

score = seu@assays$METABOLISM$score
saveRDS(score,'ssgsea.rds')

celltype = 'Myeloid'
celltype_index = curated_meta[curated_meta$celltype_major_order %in% c(celltype),] %>% rownames()
celltype_minor = names(table(curated_meta[curated_meta$celltype_major_order %in% c(celltype),'celltype_minor']))
celltype_minor_index_list = list()

celltype_score = score[,celltype_index] %>% t() %>% as.data.frame()
celltype_score$CellID = rownames(celltype_score)

sub_meta = curated_meta[c('SampleID','SampleTimepoint','celltype_minor')]
sub_meta$CellID = rownames(sub_meta)

select_celltype = 'Macrophage-TREM2'

df = left_join(celltype_score,sub_meta, by = 'CellID') %>% dplyr::filter(.,celltype_minor == select_celltype)
df$SampleTimepoint = if_else(df$SampleTimepoint == 'S1','pre','post') %>% factor(., levels = c('pre','post'))

limmstest = function(df,SampleTimepoint){
  require(limma)
  group <- df[[SampleTimepoint]] %>% as.factor()
  desigN <- model.matrix(~ 0 + group) 
  colnames(desigN) <- levels(group)
  
  fit = lmFit(df[,c(1:84)] %>% t(), desigN)
  fit2 <- eBayes(fit)
  diff=topTable(fit2,coef=2,number=Inf,adjust = 'fdr',p.value = 0.001)
}

diff = limmstest(df,'SampleTimepoint')

meta_interest = c('Fatty acid degradation')

ggplot(data = df,aes(x = SampleTimepoint, y = .data[[meta_interest]]))+
  geom_violin(aes(fill = SampleTimepoint),width = 0.7)+
  geom_boxplot(aes(fill = SampleTimepoint),width = 0.2, position = position_dodge(0.7),outlier.shape = NA)+
  scale_fill_manual(values = timepoint_color_p)+
  labs(x = '',y = 'Score', title = meta_interest)+
  stat_compare_means(comparisons = list(c('pre','post')),
                     method = "wilcox.test",
                     label = "p.format",tip.length = 0,size = 2,lwd = 0.5)+
  theme_classic()+NoLegend()+
  theme(plot.title = element_text(hjust = 0.5))

