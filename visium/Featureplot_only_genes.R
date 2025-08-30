# only feature -----------
my_feature <- c('SDC4','SDC1','CD44','ITGA11','ITGB1','ITGA1','ITGA2','ITGA3','ITGAV','ITGB8','COL4A4')
my_feature <- 'COL1A1'
for (i in 1:length(my_feature)){
  # my_number <- c(2.5, 2.5, 2, 1, 2, 2, 2, 2, 2, 2)
  # my_lim <- c(0 , my_number[i])
  # i = 1
  g.t1 <- SpatialFeaturePlot(lst[[1]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('Pri1(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  g.t2 <- SpatialFeaturePlot(lst[[2]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('Pri2(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  g.l1 <- SpatialFeaturePlot(lst[[3]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('LN1(-)')  #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  g.l2 <- SpatialFeaturePlot(lst[[4]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('LN2(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  g.t3 <- SpatialFeaturePlot(lst[[5]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('Pri3(+)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  g.l3 <- SpatialFeaturePlot(lst[[6]], feature = my_feature[i])+
    theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
    ggtitle('LN3(+)')  #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
  
  my.g <- list(g.t1, g.t2, g.l1, g.l2, g.t3, g.l3)
  assign(paste0(my_feature[i]), my.g)
}

my_visium <- list(SDC4, SDC1, CD44, ITGA11, ITGB1, ITGA1, ITGA2, ITGA3, ITGAV, ITGB8, COL4A4)
my_visium<-list(COL1A1)
save_path <- ('/home/sohee/analysis/data/HNSCC/results/seruat_visium/decided_visium_three/COLLAGEN/')

for (i in 1:length(my_feature)){
  # i = 1
  tiff(filename=paste0(save_path, paste0("Visium_",my_feature[i],".tiff")), width = 450, height = 60, unit = "mm", bg = "transparent", res = 300)
  p <- ggarrange(my_visium[[i]][[1]],my_visium[[i]][[3]], my_visium[[i]][[2]], my_visium[[i]][[4]], my_visium[[i]][[5]], my_visium[[i]][[6]], nrow = 1, ncol = 6)
  print(p)
  print(i)
  dev.off()
}
