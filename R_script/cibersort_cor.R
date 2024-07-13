library(ggpmisc)
TCGA.Kallisto.fullIDs.cibersort.relative <- read.delim("~/Project/TCGAPanCan/TCGA.Kallisto.fullIDs.cibersort.relative.tsv")
df = dplyr::filter(TCGA.Kallisto.fullIDs.cibersort.relative, CancerType == 'BRCA')
df = df[!duplicated(df$SampleID),]
rownames(df) = df$SampleID
df$SampleID = NULL
df$CancerType = NULL

cor.mat = cor(df[,c(1:22)])
cor_res = cor.mat[,'Macrophages.M2']

df[df ==0] = NA

p1 = ggplot(df,aes(x = Macrophages_M2_CIBERSORT, y = T_cells_CD8_CIBERSORT))+
  geom_point(color = 'skyblue',size = 1,na.rm = T)+
  geom_smooth(color = 'darkblue',method = 'lm',se = T)+ 
  theme_classic()+
  stat_correlation(formula = y~x, aes(label = paste(after_stat(cor.label),
                                                    after_stat(p.value.label),sep = '~~~')),
                   label.x = 0.5)+
  labs(x = 'TREM2-Mac', y = 'CD8-T-active')
                   

p2 = ggplot(df,aes(x = Macrophages_M2_CIBERSORT, y = T_cells_follicular_helper_CIBERSORT))+
  geom_point(color = 'skyblue',size = 1, na.rm = T)+
  geom_smooth(color = 'darkblue',method = 'lm',se = T,na.rm = T)+ 
  theme_classic()+
  stat_correlation(formula = y~x, aes(label = paste(after_stat(cor.label),
                                                    after_stat(p.value.label),sep = '~~~')),
                   label.x = 0.5)+
  labs(x = 'TREM2-Mac', y = 'Tfh')

p3 = ggplot(df,aes(x = Macrophages_M2_CIBERSORT, y = NK_cells_activated_CIBERSORT))+
  geom_point(color = 'skyblue',size = 1, na.rm = T)+
  geom_smooth(color = 'darkblue',method = 'lm',se = T,na.rm = T)+ 
  theme_classic()+
  stat_correlation(formula = y~x, aes(label = paste(after_stat(cor.label),
                                                    after_stat(p.value.label),sep = '~~~')),
                   label.x = 0.5)+
  labs(x = 'TREM2-Mac', y = 'NK-active')

p4 = ggplot(df,aes(x = Macrophages_M2_CIBERSORT, y = Monocytes_CIBERSORT))+
  geom_point(color = 'skyblue',size = 1, na.rm = T)+
  geom_smooth(color = 'darkblue',method = 'lm',se = T,na.rm = T)+ 
  theme_classic()+
  stat_correlation(formula = y~x, aes(label = paste(after_stat(cor.label),
                                                    after_stat(p.value.label),sep = '~~~')),
                   label.x = 0.5)+
  labs(x = 'TREM2-Mac', y = 'Treg')



p = (p1|p2)/(p3|p4)
p
ggsave('fuscc_cor_trem2.pdf',p,width = 5.69,height = 4.27,dpi = 300)
