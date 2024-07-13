
source("~/Project/MultiOmics/code/func/rungsea.R")

# enrichment res

gene_select = c('GOBP_RIBOSOME_ASSEMBLY','GOBP_TRANSLATIONAL_ELONGATION','GOBP_MRNA_TRANSCRIPTION',
                'GOBP_NUCLEAR_ENVELOPE_ORGANIZATION','GOBP_RESPONSE_TO_EPIDERMAL_GROWTH_FACTOR',
                'GOBP_RESPONSE_TO_PHEROMONE','GOBP_KERATINIZATION','GOBP_IMMUNOGLOBULIN_PRODUCTION','GOBP_REGULATION_OF_LACTATION')

plotdf = df %>% dplyr::filter(., NES>3 | NES< (-1.5) ) %>% dplyr::filter(., p.adjust<0.05)

plotdf = df[gene_select,]
plotdf$group = if_else(plotdf$NES>0,'up','down')

#plotdf$Description = gsub('_..','',plotdf$Description)
p = 
  ggplot(plotdf)+
  geom_col(aes(reorder(Description,NES),y = NES,fill=pvalue))+
  scale_fill_gsea(alpha = 0.6,reverse = T)+
  #scale_fill_manual(values = c("#8D4873","#1084A4"))+
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
  NoLegend()+
  geom_segment(aes(y=0, yend=0,x=0,xend=nrow(plotdf)+1))+
  geom_text(data = plotdf[plotdf$NES> 0.5,],aes(x=Description, y=-0.05, label=Description),
            hjust=1, size=4)+
  geom_text(data = plotdf[plotdf$NES< -0.5,],aes(x=Description, y=0.05, label=Description),
            hjust=0, size=4)+
  scale_x_discrete(expand = expansion(mult = c(0,0)))+
  scale_y_continuous(breaks = c(-1, -0.5, 0, .5, 1))+
  labs(title = '',x='', y='Normalized Enrichment Score')
p
ggsave('tf_gsea_bar.png',width = 3.18,height = 2.52,dpi = 300)
