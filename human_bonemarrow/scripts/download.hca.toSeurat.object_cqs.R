library("HCAData")
library(HCAData)
sce.bone <- HCAData('ica_bone_marrow', as.sparse=F)
sce.bone$Donor <- sub("_.*", "", sce.bone$Barcode)

library(Seurat)
sce <- CreateSeuratObject(as(assay(sce.bone, "counts"), "dgCMatrix"))
saveRDS(sce, file = "seurat_objects/hca.scedgCM.rds")




