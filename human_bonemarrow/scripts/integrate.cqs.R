#cd /rsrch3/home/lym_myl_rsch/MCL_Lab/Qingsong/PDX_single_CQS/SC145_2023/azimuth-references-master/human_bonemarrow
#conda activate R403


library(Seurat)
library(future)
library(future.apply)
#plan("multiprocess", workers = 8)
options(future.globals.maxSize = Inf)
set.seed(1234)
library(stringr)

process_objects <- function(x) {
  x$percent.mt <- PercentageFeatureSet(object = x, pattern = "^MT-")
  x <- x[, x$percent.mt < 5]
  x <- SCTransform(x, variable.features.n = 3000)
  x <- RunPCA(x)
  return(x)
}

hca <- readRDS("seurat_objects/hca.scedgCM.rds")
hca@meta.data$donor <-  str_split_fixed(rownames(hca@meta.data), "_", 2)[,1]

hca <- hca[, hca$nCount_RNA > 500]
hca <- SplitObject(object = hca, split.by = "donor")

hca <- lapply(hca, function(x) {
  x$orig.ident <- x$donor
  x
})


args <-   c( "seurat_objects/nhlbi.rds",
             "seurat_objects/hca_list.rds",
             "seurat_objects/mpal.rds",
           "seurat_objects/integrated_sct.rds")

nhlbi <- readRDS(args[1])
mpal <- readRDS(args[3])

nhlbi <- lapply(nhlbi, function(x) {
  x$orig.ident <- x$dataset
  x
})


# sct integrate
hca <- future_lapply(hca, process_objects)


nhlbi <- future_lapply(nhlbi, process_objects)

save(hca, nhlbi, file = "tmp.rds")

mpal <- UpdateSeuratObject(mpal)

DefaultAssay(mpal) <- "RNA"
obj.list <- SplitObject(mpal, split.by = "tissue")


process_objects2 <- function(x) {
  x$percent.mt <- PercentageFeatureSet(object = x, pattern = "^MT-")
  #x <- x[, x$percent.mt < 5]
  x <- SCTransform(x, variable.features.n = 3000)
  x <- RunPCA(x)
  return(x)
}
obj.list <- future_lapply(obj.list, process_objects2)


obj.list <- c(obj.list, nhlbi, hca)



ifnb.ori <- RunUMAP(obj.list)
, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")

#########################################
ifnb.ori <- NormalizeData(obj.list)
ifnb.ori <- FindVariableFeatures(ifnb.ori)
ifnb.ori <- ScaleData(ifnb.ori)
ifnb.ori <- RunPCA(ifnb.ori)

ifnb.ori <- FindNeighbors(ifnb.ori, dims = 1:30, reduction = "pca")
ifnb.ori <- FindClusters(ifnb.ori, resolution = 2, cluster.name = "unintegrated_clusters")

ifnb.ori <- RunUMAP(ifnb.ori, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")















features <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 3000)
obj.list <- PrepSCTIntegration(obj.list, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  reduction = "rpca",
  normalization.method = "SCT",
  anchor.features = features,
  dims = 1:30
)

integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)

saveRDS(object = integrated, file = args[4])

