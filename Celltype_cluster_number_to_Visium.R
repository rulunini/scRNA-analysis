# < Needed library > ----------
library(Seurat)


# < Prepare your data > ----------

my.path <- 'write_your_path'

# This is your visium dataset 
load(paste0(my.path, 'spatial_eachsample_seperated.RData'))
# This is .csv file that already created to identify enriched celltypes from divided clusters at visium 
my.df <- read.csv(paste0(my.path, 'enrichment_signif_celltype.csv'))


# < Create variables > -----------

# Convert individual visium data to list format
my.visium.lst <- list(tm1, ln1, tm2, ln2, tm3, ln3)

my.celltype <- colnames(my.df)
my.celltype <- paste0(my.celltype, '.enrich')


# < Just run the code! > ----------

for (i in 1:length(my.visium.lst)){
  
  for (z in 1:length(my.celltype)){

  # Create new metadata that include my resolution
  my.visium.lst[[i]]$new <- my.visium.lst[[i]]$integrated_snn_res.1
  my.visium.lst[[i]]$new <- my.visium.lst[[i]]$integrated_snn_res.1[my.visium.lst[[i]]$integrated_snn_res.1 %in% as.numeric(na.omit(my.df[,z]))]
  colnames(my.visium.lst[[i]]@meta.data)[length(colnames(my.visium.lst[[i]]@meta.data))] <- my.celltype[z]
  
  
  my.visium.lst[[i]]@meta.data[[my.celltype[z]]][my.visium.lst[[i]]@meta.data[[my.celltype[z]]] %in% as.numeric(na.omit(my.df[,z]))] <- '1'
  my.visium.lst[[i]]@meta.data[[my.celltype[z]]][my.visium.lst[[i]]@meta.data[[my.celltype[z]]] %in% NA] <- '0'
  my.visium.lst[[i]]@meta.data[[my.celltype[z]]] <- factor(my.visium.lst[[i]]@meta.data[[my.celltype[z]]])
  Idents(my.visium.lst[[i]]) <- my.celltype[z]
  
  }
}




# < Save to RData > ----------

save.path <- 'write_your_path'
save(my.visium.lst, file = paste0(save.path, '/spatial_sixsample_integration_celltype.RData'))

