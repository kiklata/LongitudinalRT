
stackplot = function(obj, ylabuse = FALSE, mycol= c("#8DD3C7", "#FFFFB3", "#80B1D3", "#8DA0CB", "#B3DE69",
                                                    "#FB8072", "#BEBADA", "#D9D9D9", "#FC8D62", "#CCEBC5",
                                                    "#FCCDE5")){
  p <- ggplot(data = obj) +
    geom_col(aes(x = xname,y = value/CD8,fill = variable),
             position = position_stack(),#stat = 'identity',
             width = 0.9) +
    scale_fill_manual(values = mycol) +
    #scale_y_continuous(labels = scales::percent_format()) +
    theme_gray(base_size = 18) + guides(fill = guide_legend(title = NULL,label.position = 'right',nrow = 2))+
    xlab('')
  if (ylabuse == T) {
    p = p + ylab('Proportion') 
  }
  if (ylabuse == F) {
    p = p + ylab('') + theme(axis.ticks.y = element_blank())
  }
  p = p + 
    theme(axis.text.x = element_text(family = 'arial',size = 10,angle = 45,vjust = 1,hjust = 0.5,colour = 'black'),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_blank(),
          legend.position = 'bottom',legend.direction = 'horizontal'
          #plot.margin = margin(t = 3,b = 3,r = 1,unit = 'cm')
    )
  return(p)
}

anno = read.delim('~/Project/MultiOmics/data/skin/res/cellbender_celltype.csv',row.names = 1,sep = ',')
anno = new_anno@meta.data
anno$cl_comp = if_else(anno$cl_major %in% c('CD4-T','CD8-T','NKT'),'T cell',
                       if_else(anno$cl_major %in% c('Macrophage','Monocyte','DC'),'Myeloid',
                               if_else(anno$cl_major %in% c('Pericyte','vSMC'),'PVL',
                                       if_else(anno$cl_major %in% c('iCAF','mCAF','vCAF'),'CAF',
                                               anno$cl_major))))  
saveRDS(anno, file = 'anno.rds')
#anno = dplyr::filter(anno, cl_major == 'Keratinocyte')

ptexpan = as.data.frame(table(anno$pseudo_sample,anno$cl_comp))
ptexpan = tidyr::spread(ptexpan,key = 'Var2',value = 'Freq')
ptexpan$total = apply(ptexpan[,2:ncol(ptexpan)],1,sum)

celltype_list = c("KC-Basal-COL17A1","KC-Basal-ITGA6", "KC-Suprabasal-DSC3" ,    
                  "KC-Suprabasal-KRT6A","KC-Suprabasal-LYPD3","KC-Suprabasal-PLD1",
                  "KC-Spinous-AZGP1","KC-Spinous-SPINK5",           
                  "KC-Proliferating-DIAPH3","KC-Cycling" )

stackbar_sp_order = c('P01_pre','P02_pre','P03_pre','P01_post','P02_post','P03_post')

stackbar_cl_major_order = c('CD4-T','CD8-T','NKT','B cell','Macrophage','Monocyte','DC','Mast',
                      'Endothelial','Pericyte','vSMC',
                      'Epithelial',
                      'iCAF','mCAF','vCAF')
stackbar_cl_comp_order = c('T cell','B cell','Myeloid','Mast','Endothelial','PVL','Epithelial','CAF')

df = ptexpan %>% reshape2::melt(.,c('Var1','total'))

df$variable = factor(df$variable,levels = stackbar_cl_comp_order)
df$Var1 = factor(df$Var1, levels = rev(stackbar_sp_order))

p = ggplot(data = df) +
  geom_col(aes(x = Var1,y = value/total,fill = variable),
           position = position_stack(),
           width = 0.9,color = 'grey25') + 
  coord_flip()+ # cl on, sample off
  scale_fill_manual(values = compart_color_p) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 6,angle = 0,hjust = 1,colour = 'black'),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 8),
        axis.text.y = element_text(size = 6, colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2,'cm'))+ # cl 0.2, sample 0.5
  labs(x = '', y = 'Proportion')+guides(fill = guide_legend(ncol= 1))
p
ggsave('sample_stackbar.pdf',p,width = 6.17,height = 1.83,dpi = 300)
ggsave('cl_minor_stackbar.png',p,width = 5.66,height = 1.71,dpi = 300)
