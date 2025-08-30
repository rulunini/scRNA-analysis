# < Needed library > ----------
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(ggpubr)

# < Prepare your data > ----------

my.path <- 'write_your_path'
load(paste0(my.path, 'spatial_sixsample_integration_celltype.RData'))


# < Visualize as tumor > ----------

my.color <- c('#EAEAEA','#E82D23')
for (i in 1:length(my.visium.lst[[i]])){
  Idents(my.visium.lst[[i]]) <- 'write_the_celltypes_that_you_want'
}

a1 <- SpatialDimPlot(my.visium.lst[[1]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri1(-)')
a2 <- SpatialDimPlot(my.visium.lst[[2]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN1(-)')
a3 <- SpatialDimPlot(my.visium.lst[[3]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a4 <- SpatialDimPlot(my.visium.lst[[4]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN2(-)')
a5 <- SpatialDimPlot(my.visium.lst[[5]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri3(+)')
a6 <- SpatialDimPlot(my.visium.lst[[6]], cols = my.color)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN3(+)')

ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)


# < Save .tiff > ----------

save.path = 'write_your_path'
tiff(filename = paste0(my.path, "MIA_spatialdimplot_tumor.tiff"), 
     width = 330, height = 60, unit = "mm", bg = "transparent", res = 300)
ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)
dev.off()