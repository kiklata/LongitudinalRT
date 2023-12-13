#devtools::install_github("xl0418/ggradar2",dependencies=TRUE)
library(ggradar2)

obj = CD8.Tex.harmony@meta.data

plotdf = obj[,c(26:30,32)]
plotdf = aggregate(plotdf[,c(1:5)],FUN = mean,by = list(plotdf[,'manual.celltype.Tex']))

colnames(plotdf)[1] = 'group'

plotdf[2:6] = apply(
  plotdf[2:6],2,
  FUN = function(data) {
    new.data = (data - min(data)) / (max(data) - min(data))
  }
)
mycolor = c('#4dbbd5','#3c5488','#f39b7f','#00a087')

ggradar2(plotdf, plot.legend = T, grid.label.size = F, 
         group.colours = mycolor)

