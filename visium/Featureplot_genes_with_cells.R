
# gene * cellcluster 
for (i in 1:6){
  
  # TNFSF4-TNFRSF4 -> OX40
  lst[[i]]$TNFSF4.ligand <- lst[[i]]@assays$SCT@data['TNFSF4',] * 
    (as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$ImmatureT.enrich)) +
       as.numeric(as.character(lst[[i]]$Endo.enrich)))
  Idents(lst[[i]])<-'TNFSF4.ligand'
  
  lst[[i]]$TNFRSF4.receptor <- lst[[i]]@assays$SCT@data['TNFRSF4',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich))
  Idents(lst[[i]])<-'TNFRSF4.receptor'
  
  # JAG1-NOTCH21 -> NOTCH
  lst[[i]]$JAG1.ligand <- lst[[i]]@assays$SCT@data['JAG1',] * 
    (as.numeric(as.character(lst[[i]]$Endo.enrich))+
       as.numeric(as.character(lst[[i]]$ImmatureT.enrich)) +
       as.numeric(as.character(lst[[i]]$Mast.enrich)) +
       as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$Myo.enrich)) +
       as.numeric(as.character(lst[[i]]$NAF.enrich)) +
       as.numeric(as.character(lst[[i]]$Normal.enrich)))
  Idents(lst[[i]])<-'JAG1.ligand'
  # JAG2
  lst[[i]]$JAG2.ligand <- lst[[i]]@assays$SCT@data['JAG2',] * 
    (as.numeric(as.character(lst[[i]]$Endo.enrich))+
       as.numeric(as.character(lst[[i]]$ImmatureT.enrich)) +
       as.numeric(as.character(lst[[i]]$Mast.enrich)) +
       as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$Myo.enrich)) +
       as.numeric(as.character(lst[[i]]$NAF.enrich)) +
       as.numeric(as.character(lst[[i]]$Normal.enrich)))
  Idents(lst[[i]])<-'JAG2.ligand'
  lst[[i]]$NOTCH2.receptor <- lst[[i]]@assays$SCT@data['NOTCH2',] * 
    (as.numeric(as.character(lst[[i]]$Endo.enrich))+
       as.numeric(as.character(lst[[i]]$ImmatureT.enrich)) +
       as.numeric(as.character(lst[[i]]$Mast.enrich)) +
       as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$Myo.enrich)) +
       as.numeric(as.character(lst[[i]]$NAF.enrich)) +
       as.numeric(as.character(lst[[i]]$Normal.enrich)))
  Idents(lst[[i]])<-'NOTCH2.receptor'
  
  # CSF3-CSF3R -> CSF
  lst[[i]]$CSF3.ligand <- lst[[i]]@assays$SCT@data['CSF3',] * 
    (as.numeric(as.character(lst[[i]]$Endo.enrich))+
       as.numeric(as.character(lst[[i]]$ImmatureT.enrich)) +
       as.numeric(as.character(lst[[i]]$Mast.enrich)) +
       as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$Tumor.enrich)))
  Idents(lst[[i]])<-'CSF3.ligand'
  lst[[i]]$CSF3R.receptor <- lst[[i]]@assays$SCT@data['CSF3R',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich))
  # as.numeric(as.character(lst[[i]]$Peicyte.enrich)) *
  Idents(lst[[i]])<-'CSF3R.receptor'
  
  # EFNA1-EPHA2, EFNA3-EPHA2-> EPHA
  lst[[i]]$EFNA1.ligand <- lst[[i]]@assays$SCT@data['EFNA1',] * 
    (as.numeric(as.character(lst[[i]]$NAF.enrich)) +
       as.numeric(as.character(lst[[i]]$Myo.enrich)))
  Idents(lst[[i]])<-'EFNA1.ligand'
  lst[[i]]$EFNA3.ligand <- lst[[i]]@assays$SCT@data['EFNA3',] * 
    (as.numeric(as.character(lst[[i]]$NAF.enrich)) +
       as.numeric(as.character(lst[[i]]$Myo.enrich)))
  Idents(lst[[i]])<-'EFNA3.ligand'
  lst[[i]]$EPHA2.receptor <- lst[[i]]@assays$SCT@data['EPHA2',] * 
    as.numeric(as.character(lst[[i]]$Macro.enrich))
  Idents(lst[[i]])<-'EPHA2.receptor'
  
  # SELE-GLG1 -> SELPLG
  lst[[i]]$SELE.ligand <- lst[[i]]@assays$SCT@data['SELE',] * 
    (as.numeric(as.character(lst[[i]]$DC.enrich))+
       as.numeric(as.character(lst[[i]]$CAF.enrich)) )
  # as.numeric(as.character(lst[[i]]$Pericyte.enrich)) *
  Idents(lst[[i]])<-'SELE.ligand'
  lst[[i]]$GLG1.receptor <- lst[[i]]@assays$SCT@data['GLG1',] * 
    (as.numeric(as.character(lst[[i]]$CAF.enrich))+
       as.numeric(as.character(lst[[i]]$DC.enrich)) +
       as.numeric(as.character(lst[[i]]$Macro.enrich))+
       as.numeric(as.character(lst[[i]]$Myo.enrich)))
  # as.numeric(as.character(lst[[i]]$Pericyte.enrich))
  Idents(lst[[i]])<-'GLG1.receptor'
}

# TNFSF4-TNFRSF4 -> OX40
# JAG1-NOTCH21 -> NOTCH
# CSF3-CSF3R -> CSF
# EFNA1-EPHA2, EFNA3-EPHA2-> EPHA
# SELE-GLG1 -> SELPLG

# TNF -> TNFRSF1B, TNFRSF1A
for (i in 1:6){
  lst[[i]]$TNF.CAF <- lst[[i]]@assays$SCT@data['TNF',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich)) 
  Idents(lst[[i]])<-'TNF.CAF'
  lst[[i]]$TNF.Tumor <- lst[[i]]@assays$SCT@data['TNF',] * 
    as.numeric(as.character(lst[[i]]$Tumor.enrich)) 
  Idents(lst[[i]])<-'Tumor'
}

# TWEAK -> TNFSF12, TNFRSF12A
for (i in 1:6){
  # TNF -> TNFSF12, TNFRSF12A
  lst[[i]]$TNFSF12.CAF <- lst[[i]]@assays$SCT@data['TNFSF12',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich)) 
  Idents(lst[[i]])<-'TNFSF12.CAF'
  lst[[i]]$TNFSF12.Pericyte <- lst[[i]]@assays$SCT@data['TNFSF12',] * 
    as.numeric(as.character(lst[[i]]$Pericyte.enrich)) 
  Idents(lst[[i]])<-'TNFSF12.Peicyte'
}
# TGFB1-ACVR1+ACVR1B+TGFBR11+TGFBR21
for (i in 1:6){
  # TGFB1-ACVR1+ACVR1B+TGFBR11+TGFBR21
  lst[[i]]$TGFB1.Tumor <- lst[[i]]@assays$SCT@data['TGFB1',] * 
    as.numeric(as.character(lst[[i]]$Tumor.enrich)) 
  Idents(lst[[i]])<-'TGFB1.Tumor'
  
  lst[[i]]$TGFB1.Mast <- lst[[i]]@assays$SCT@data['TGFB1',] * 
    as.numeric(as.character(lst[[i]]$Mast.enrich)) 
  Idents(lst[[i]])<-'TGFB1.Mast'
  
  lst[[i]]$TGFB1.NAF <- lst[[i]]@assays$SCT@data['TGFB1',] * 
    as.numeric(as.character(lst[[i]]$NAF.enrich)) 
  Idents(lst[[i]])<-'TGFB1.NAF'
}  

# EPHA -> EFNA1, EPHA2, EPHA3
for (i in 1:6){
  lst[[i]]$EPHA1.NAF <- lst[[i]]@assays$SCT@data['EFNA1',] * 
    as.numeric(as.character(lst[[i]]$NAF.enrich)) 
  Idents(lst[[i]])<-'EFNA1.NAF'
  lst[[i]]$EPHA1.Myo <- lst[[i]]@assays$SCT@data['EFNA1',] * 
    as.numeric(as.character(lst[[i]]$Myo.enrich)) 
  Idents(lst[[i]])<-'EFNA1.Myo'
}

# CDH1 -> CDH1, ITGA1, ITGB1
for (i in 1:6){
  lst[[i]]$CDH1.NAF <- lst[[i]]@assays$SCT@data['CDH1',] * 
    as.numeric(as.character(lst[[i]]$NAF.enrich)) 
  Idents(lst[[i]])<-'CDH1.NAF'
  lst[[i]]$ITGA1.CAF <- lst[[i]]@assays$SCT@data['ITGA1',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich)) 
  Idents(lst[[i]])<-'ITGA1.CAF'
  lst[[i]]$ITGB1.CAF <- lst[[i]]@assays$SCT@data['ITGB1',] * 
    as.numeric(as.character(lst[[i]]$CAF.enrich)) 
  Idents(lst[[i]])<-'ITGB1.CAF'
}

# SELE -> CDH1, ITGA1, ITGB1
for (i in 1:6){
  lst[[i]]$SELE.Myo <- lst[[i]]@assays$SCT@data['SELE',] * 
    as.numeric(as.character(lst[[i]]$Myo.enrich)) 
  Idents(lst[[i]])<-'SELE.Myo'
}
my_feature <- c('TNFSF4.ligand','TNFRSF4.receptor', 'JAG1.ligand','JAG2.ligand','NOTCH2.receptor')
my_feature <- c('CSF3.ligand','CSF3R.receptor','EFNA1.ligand','EFNA3.ligand','EPHA2.receptor','SELE.ligand','GLG1.receptor')

# for (i in 1:length(my_feature)){
my_number <- c(1, 1, 3, 3, 1)
# my_lim <- c(0 , my_number[i])

g.t1 <- SpatialFeaturePlot(lst[[1]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Pri1(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
g.t2 <- SpatialFeaturePlot(lst[[2]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Pri2(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
g.l1 <- SpatialFeaturePlot(lst[[3]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('LN1(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
g.l2 <- SpatialFeaturePlot(lst[[4]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('LN2(-)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
g.t3 <- SpatialFeaturePlot(lst[[5]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('Pri3(+)')# +scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))
g.l3 <- SpatialFeaturePlot(lst[[6]], feature = my_feature[i])+
  theme(legend.position = 'right', legend.title = element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle('LN3(+)') #+scale_fill_gradientn(limits = my_lim, breaks = c(0,max(my_lim),max(my_lim)/2),colors =c('#DBDBDB','#E2B126','#FF0909','#CC0000'))

my.g <- list(g.t1, g.t2, g.l1, g.l2, g.t3, g.l3)
assign(paste0(my_feature[i]), my.g)
}

my_visium <- list(TNFSF4.ligand, TNFRSF4.receptor, JAG1.ligand, JAG2.ligand, NOTCH2.receptor)
my_name <- c('TNFSF4_ligand', 'TNFRSF4_receptor','JAG1_ligand','JAG2_ligand','NOTCH2_receptor')

my_visium <- list(CSF3.ligand, CSF3R.receptor, EFNA1.ligand, EFNA3.ligand, EPHA2.receptor, SELE.ligand, GLG1.receptor)
my_name <- c('CSF3_ligand', 'CSF3R_receptor','EFNA1_ligand','EFNA3_ligand','EPHA2_receptor','SELE_ligand','GLG1_receptor')

i = 5
tiff(filename=paste0(save_path, paste0("Visium_",my_name[i],".tiff")), width = 450, height = 60, unit = "mm", bg = "transparent", res = 300)
ggarrange(my_visium[[i]][[1]],my_visium[[i]][[3]], my_visium[[i]][[2]], my_visium[[i]][[4]], my_visium[[i]][[5]], my_visium[[i]][[6]], nrow = 1, ncol = 6)
dev.off()

