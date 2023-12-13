
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
