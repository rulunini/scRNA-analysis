# < Needed library > ----------
library(Seurat)
library(SeuratObject)


# < Prepare your data > ----------

load('/home/sohee/analysis/data/HNSCC/mydata/mia/spatial_eachsample_seperated.RData')
load('/home/sohee/analysis/data/HNSCC/mydata/mia/spatial_eightsample_integration.RData')


# < Visualize as tumor > ----------
visium.int$tumor_cluster <- visium.int$integrated_snn_res.1[visium.int$integrated_snn_res.1%in% c(0,2,4,9,11,12,13)]
visium.int$tumor_cluster <- factor(visium.int$tumor_cluster)
visium.int$tumor_cluster <- as.character(visium.int$tumor_cluster)
visium.int$tumor_cluster[visium.int$tumor_cluster%in%c(0,2,4,9,11,12,13)] <- 'tumor'
visium.int$tumor_cluster[visium.int$tumor_cluster%in%NA] <- 'others'
visium.int$tumor_cluster <- factor(visium.int$tumor_cluster)


color = c('#EAEAEA','#E82D23')
Idents(visium.int) <- 'tumor_cluster'
a1 <- SpatialDimPlot(visium.int, images = 'TM1', cols = color, pt.size.factor = 4+2.2-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri1(-)')
a2 <- SpatialDimPlot(visium.int, images = 'LN1', cols = color, pt.size.factor = 4+1.9-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN1(-)')
a3 <- SpatialDimPlot(visium.int, images = 'TM2', cols = color, pt.size.factor = 4+2.3-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a4 <- SpatialDimPlot(visium.int, images = 'LN2', cols = color, pt.size.factor = 4+2.1-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN2(-)')
a5 <- SpatialDimPlot(visium.int, images = 'TM3', cols = color, pt.size.factor = 4+2.2-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri3(+)')
a6 <- SpatialDimPlot(visium.int, images = 'LN3', cols = color, pt.size.factor = 4+1.9-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN3(+)')
ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)


# < Save .tiff > ----------

tiff(filename=paste0(save_path, "MIA_spatialdimplot_tumor.tiff"), 
     width = 330, height = 60, unit = "mm", bg = "transparent", res = 300)
ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)
dev.off()