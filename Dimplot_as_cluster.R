# < Needed library > ----------
library(Seurat)


# < Prepare your data > ----------

my.path <- '/home/sohee/analysis/data/HNSCC/mydata/mia/' #'write_your_path'
load(paste0(my.path, 'spatial_eightsample_integration.RData'))


# < Visualize as barplot > ----------

visium.int$site2 <- visium.int$site
visium.int$site2 <- factor(visium.int$site2)
levels(visium.int$site2)[levels(visium.int$site2)=='Primary']<-'Pri'
levels(visium.int$site2)[levels(visium.int$site2)=='LymphNode']<-'LN'

visium.int$HPV2 <- visium.int$HPV
visium.int$HPV2 <- factor(visium.int$HPV2)
levels(visium.int$HPV2)[levels(visium.int$HPV2)=='Negative']<-'(-)'
levels(visium.int$HPV2)[levels(visium.int$HPV2)=='Positive']<-'(+)'

visium.int$sample2 <- visium.int$sample
visium.int$sample2 <- factor(visium.int$sample2)
levels(visium.int$sample2)[levels(visium.int$sample2)=='PT1']<-'1'
levels(visium.int$sample2)[levels(visium.int$sample2)=='PT2']<-'2'
levels(visium.int$sample2)[levels(visium.int$sample2)=='PT3']<-'3'

visium.int$label <- paste0(visium.int$site2, visium.int$sample2, visium.int$HPV2)
visium.int$label <- factor(visium.int$label)
visium.int$label <- factor(visium.int$label, levels = c('LN3(+)','Pri3(+)','LN2(-)','Pri2(-)','LN1(-)', 'Pri1(-)'))

table(visium.int$integrated_snn_res.1)
table(visium.int$label)

df <- round(prop.table(table(visium.int$label, visium.int$integrated_snn_res.1),2)*100,2)
df <- df[,colnames(df) %in% c(4, 5, 12, 15, 3, 8, 9 )]
df <- melt(df)
colnames(df) <- c('label','cluster','prop')
df$cluster <- factor(df$cluster)
df$cluster <- factor(df$cluster, levels = c(4,5, 12, 15, 3, 8, 9))
levels(df$cluster)[levels(df$cluster)==4]<-'cluster4'
levels(df$cluster)[levels(df$cluster)==5]<-'cluster5'
levels(df$cluster)[levels(df$cluster)==12]<-'cluster12'
levels(df$cluster)[levels(df$cluster)==15]<-'cluster15'
levels(df$cluster)[levels(df$cluster)==3]<-'cluster3'
levels(df$cluster)[levels(df$cluster)==8]<-'cluster8'
levels(df$cluster)[levels(df$cluster)==9]<-'cluster9'
df <- data.frame(df)
df$enrichment <- NA
df[df$cluster%in% c('cluster4','cluster5','cluster12','cluster15'),'enrichment'] <- 'HPV(-)up'
df[df$cluster%in% c('cluster3','cluster8','cluster9'),'enrichment'] <- 'HPV(+)up'

a <- df[df$enrichment == 'HPV(-)up',]
b <- df[df$enrichment == 'HPV(+)up',]
X_color <- c("#EDF151", "#FF725C","#70E0E6","#58EAAF", "#ff28a7", "#4B0082")

color <- c("#FFB400", "#FF4500","#006400","#3232ff","#ff50cf", "#C12Dff","#4B0082", 'grey')
b<-b[b$cluster %in% 'cluster8',]
g2 <-ggplot(b, aes(x = prop, y = label, fill = cluster)) +
  geom_bar(position = position_dodge(), stat = 'identity', width = 0.8, col = 'black', alpha = 0.7)+
  theme_bw()+
  facet_grid(~cluster, scales = 'free')+
  theme(legend.position = "",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(color = 'black'),
        axis.text.y = element_text(color = 'black'),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5))+
  ggtitle('')+
  scale_fill_manual(values = color[9])
g2

tiff(filename=paste0(save_path, "MIA_proportion_neg.tiff"), width = 155, height = 45, unit = "mm", bg = "transparent", res = 300)
g1
dev.off()
tiff(filename=paste0(save_path, "MIA_proportion_pos.tiff"), width = 83, height = 45, unit = "mm", bg = "transparent", res = 300)
g2
dev.off()


# Visualize as hpv spatial----------------------------------------------------------------------------------------------------
# "#E8BA5A", "#FF725C","#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA","#E82D23"

Idents(visium.int) <- 'integrated_snn_res.1'
visium.int$neg_cluster <- visium.int$integrated_snn_res.1[visium.int$integrated_snn_res.1%in% c(4, 5, 12, 15)]
visium.int$neg_cluster <- factor(visium.int$neg_cluster, levels = c(4, 5, 12, 15))
visium.int$neg_cluster <- as.character(visium.int$neg_cluster)
# visium.int$neg_cluster[visium.int$neg_cluster%in%c(4, 5, 12, 15)] <- 'HPV(-)'
visium.int$neg_cluster[visium.int$neg_cluster%in%NA] <- 'others'
visium.int$neg_cluster <- factor(visium.int$neg_cluster, levels = c(4, 5, 12, 15, 'others'))
visium.int$neg_cluster <- factor(visium.int$neg_cluster)

visium.int$pos_cluster <- visium.int$integrated_snn_res.1[visium.int$integrated_snn_res.1%in% c(3, 9)]
visium.int$pos_cluster <- as.character(visium.int$pos_cluster)
# visium.int$pos_cluster[visium.int$pos_cluster%in%c(3,9)] <- 'HPV(+)'
visium.int$pos_cluster[visium.int$pos_cluster%in%NA] <- 'others'
visium.int$pos_cluster <- factor(visium.int$pos_cluster)
visium.int$pos_cluster <- factor(visium.int$pos_cluster, levels = c(3, 9, 'others'))


visium.int$hpv_cluster <- visium.int$integrated_snn_res.1[visium.int$integrated_snn_res.1%in% c(4, 5, 12, 15, 3, 9)]
visium.int$hpv_cluster <- as.character(visium.int$hpv_cluster)
visium.int$hpv_cluster[visium.int$hpv_cluster%in%c(4, 5, 12, 15)] <- 'HPV(-)'
visium.int$hpv_cluster[visium.int$hpv_cluster%in%c(3,9)] <- 'HPV(+)'
visium.int$hpv_cluster[visium.int$hpv_cluster%in%NA] <- 'others'
visium.int$hpv_cluster <- factor(visium.int$hpv_cluster)
visium.int$hpv_cluster <- factor(visium.int$hpv_cluster, levels = c('HPV(-)','HPV(+)','others'))

Idents(visium.int) <- 'neg_cluster'
color <- c("#FFB400", "#FF4500","#006400","#3232ff","#d8bfd8" )
a1 <- SpatialDimPlot(visium.int, images = 'TM1', cols = color, pt.size.factor = 2.5-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri1(-)')
a2 <- SpatialDimPlot(visium.int, images = 'LN1', cols = color, pt.size.factor = 2.2-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN1(-)')
a3 <- SpatialDimPlot(visium.int, images = 'TM2', cols = color, pt.size.factor = 2.7-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a4 <- SpatialDimPlot(visium.int, images = 'LN2', cols = color, pt.size.factor = 2.4-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN2(-)')
a5 <- SpatialDimPlot(visium.int, images = 'TM3', cols = c(color[1:3], color[5]), pt.size.factor = 2.6-0.5, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri3(+)')
a6 <- SpatialDimPlot(visium.int, images = 'LN3', cols = color, pt.size.factor = 2.2-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN3(+)')
a <- ggarrange(a1,a2, a3, a4, a5, a6, nrow = 3, ncol = 2)
a

Idents(visium.int) <- 'pos_cluster'
color <- c( "#ff50cf", "#C12Dff","#d8bfd8") #-> ÆÄ¶û, ÇÖÇÎÅ©
a1 <- SpatialDimPlot(visium.int, images = 'TM1', cols = color, pt.size.factor = 2.5-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri1(-)')
a2 <- SpatialDimPlot(visium.int, images = 'LN1', cols = color, pt.size.factor = 2.2-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN1(-)')
a3 <- SpatialDimPlot(visium.int, images = 'TM2', cols = color, pt.size.factor = 2.7-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a4 <- SpatialDimPlot(visium.int, images = 'LN2', cols = color, pt.size.factor = 2.4-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN2(-)')
a5 <- SpatialDimPlot(visium.int, images = 'TM3', cols = color, pt.size.factor = 2.6-0.5, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a6 <- SpatialDimPlot(visium.int, images = 'LN3', cols = color, pt.size.factor = 2.2-0.5, stroke = 0.2, alpha = 0.7)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN3(+)')
b <- ggarrange(a1,a2, a3, a4, a5, a6, nrow = 3, ncol = 2)
b

tiff(filename=paste0(save_path, "MIA_spatialdimplot.tiff"), width = 190, height = 160, unit = "mm", bg = "transparent", res = 300)
ggarrange(a, b, nrow = 1)
dev.off()

a1 <- SpatialDimPlot(visium.int, images = 'TM1', cols = c("#E8BA5A", "#FF725C","#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+2.2-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri1(-)')
a2 <- SpatialDimPlot(visium.int, images = 'LN1', cols = c("#E8BA5A","#FF725C", "#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+1.9-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN1(-)')
a3 <- SpatialDimPlot(visium.int, images = 'TM2', cols = c("#E8BA5A","#FF725C", "#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+2.3-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri2(-)')
a4 <- SpatialDimPlot(visium.int, images = 'LN2', cols = c("#E8BA5A","#FF725C", "#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+2.1-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN2(-)')
a5 <- SpatialDimPlot(visium.int, images = 'TM3', cols = c("#E8BA5A", "#FF725C","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+2.2-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('Pri3(+)')
a6 <- SpatialDimPlot(visium.int, images = 'LN3', cols = c("#E8BA5A","#FF725C", "#58EAAF","#FF90A4", "#D4A1FE","#EAEAEA"), pt.size.factor = 4+1.9-0.5)+
  theme(legend.position = '', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+ggtitle('LN3(+)')
ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)


# < Save .tiff > ----------

save.path <- # 'write_your_path'
  
tiff(filename = paste0(save.path, "MIA_spatialdimplot3.tiff"), 
     width = 330, height = 60, unit = "mm", bg = "transparent", res = 300)
ggarrange(a1,a2, a3, a4, a5, a6, nrow = 1, ncol = 6)
dev.off()




# save(visium.int, file = '/home/sohee/analysis/data/HNSCC/mydata/mia/spatial_eightsample_integration_v2.RData')
# load('/home/sohee/analysis/data/HNSCC/mydata/mia/spatial_eightsample_integration_v2.RData')