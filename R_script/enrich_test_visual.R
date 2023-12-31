## ssGSEA

# RUN scMetabolism
curated_meta = read.delim('~/Project/MultiOmics/data/skin/res/cellbender_celltype.csv',row.names = 1,sep = ',')
scores = rbind(h_score,meta_score) %>% t()

seu = readRDS("~/Project/MultiOmics/data/skin/res/cellbender_anno_count.rds")

score = seu@assays$METABOLISM$score
saveRDS(score,'ssgsea.rds')

celltype = 'Keratinocyte'
celltype_index = curated_meta[curated_meta$cl_major %in% c(celltype),] %>% rownames()
celltype_minor = names(table(curated_meta[curated_meta$cl_major %in% c(celltype),'cl_minor']))
celltype_minor_index_list = list()

celltype_score = scores[celltype_index,] %>% as.data.frame()
celltype_score$CellID = rownames(celltype_score)

sub_meta = curated_meta[c('SampleID','SampleType','cl_minor','cl_subset')]
sub_meta$CellID = rownames(sub_meta)

select_celltype = 'KC-Basal-COL17A1'

df = left_join(celltype_score,sub_meta, by = 'CellID') %>% dplyr::filter(., cl_subset == select_celltype)
df$SampleType = if_else(df$SampleType == 'H','ARD','Normal') %>% factor(., levels = c('ARD','Normal'))

limmstest = function(df,SampleType){
  require(limma)
  group <- df[[SampleType]] %>% as.factor()
  desigN <- model.matrix(~ 0 + group) 
  colnames(desigN) <- levels(group)
  
  fit = lmFit(df[,c(1:134)] %>% t(), desigN)
  fit2 <- eBayes(fit)
  diff=topTable(fit2,coef=2,number=Inf,adjust = 'fdr',p.value = 0.001)
}

diff = limmstest(df,'SampleType')

meta_interest = c('HALLMARK_IL6_JAK_STAT3_SIGNALING')

#library(ggpubr)
ggplot(data = df,aes(x = SampleType, y = .data[[meta_interest]]))+
  geom_violin(aes(fill = SampleType),width = 0.7)+
  geom_boxplot(aes(fill = SampleType),width = 0.2, position = position_dodge(0.7),outlier.shape = NA)+
  #scale_fill_manual(values = timepoint_color_p)+
  labs(x = '',y = 'Score', title = meta_interest)+
  stat_compare_means(comparisons = list(c('ARD','Normal')),
                     method = "wilcox.test",
                     label = "p.format",tip.length = 0,size = 2,lwd = 0.5)+
  theme_classic()+
  NoLegend()+
  theme(plot.title = element_text(hjust = 0.5))

