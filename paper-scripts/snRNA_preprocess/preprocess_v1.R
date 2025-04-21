pkgs = c("Seurat", "sctransform", "glmGamPoi",  "DoubletFinder",
         "dplyr", "plyr", "tibble", "fs", "ggplot2", "parallel", "future")
lapply(pkgs, library, character.only = TRUE)

# read_file = function(f){
#   dat.fl = readRDS(f)}

# 필터 정의
min_cells = 3 
min_features = 200 
max_percent_mt = 10 
min_nFeature_RNA = 500
max_nFeature_RNA = 10000
min_nCount_RNA = 1000
max_nCount_RNA = 30000

# 파일 경로 탐색
get_files = function(root_path) {
  files = dir_ls(root_path, recurse = TRUE, regexp = "filtered_feature_bc_matrix.h5$")
  return(files)
}

# file_path = "/Volumes/WindySSD/YUGBM_OC/cellranger/run_cellranger_count/CN8_RawSVZ_Nuclei_control1/outs/filtered_feature_bc_matrix.h5"

# 수행할 작업 정의
make_seurat_obj = function(file_path) {
  message(paste("Start with sample: ", basename(dirname(dirname(file_path)))))
  seurat_obj = CreateSeuratObject(counts = Read10X_h5(file_path, use.names = TRUE),
                                min.cells = min_cells, min.features = min_features,
                                project = "YUGBM")
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern="^MT-")
  
  p_before_filtering = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    labs(caption = "Before filtering")
  
  seurat_obj <- subset(seurat_obj, 
                       subset = 
                         percent.mt < max_percent_mt &
                         nFeature_RNA > min_nFeature_RNA & 
                         nFeature_RNA < max_nFeature_RNA &
                         nCount_RNA > min_nCount_RNA & 
                         nCount_RNA < max_nCount_RNA )
  
  p_after_filtering = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
    labs(caption = "After filtering")
  
  # SC transform
  # seurat_obj = SCTransform(seurat_obj, vars.to.regress = "percent.mt", 
  #                          return.only.var.genes = TRUE, verbose = FALSE)
  seurat_obj = SCTransform(seurat_obj, vars.to.regress = "percent.mt", 
                           return.only.var.genes = FALSE, verbose = FALSE)
  seurat_obj = RunPCA(seurat_obj, verbose = FALSE)
  
  PCpercentage = seurat_obj@reductions$pca@stdev / sum(seurat_obj@reductions$pca@stdev) * 100
  PCcumulative = cumsum(PCpercentage)
  pc1 = which(PCcumulative > 90 & PCpercentage < 5)[1]
  pc2 = sort(which((PCpercentage[1:length(PCpercentage)-1] - PCpercentage[2:length(PCpercentage)]) > 0.1),  decreasing = T)[1] + 1
  nPC = min(pc1, pc2)
  cat('selected PC value :', nPC, "\n")
  nPC = max(nPC, 10)
  
  seurat_obj = RunUMAP(seurat_obj, dims = 1:nPC, verbose = FALSE)
  seurat_obj = RunTSNE(seurat_obj, dims = 1:nPC, verbose = FALSE)
  
  seurat_obj = FindNeighbors(seurat_obj, dims = 1:nPC, verbose = FALSE)
  seurat_obj = FindClusters(seurat_obj, verbose = FALSE)

  sample = basename(dirname(dirname(file_path)))
  p_cluster_umap = DimPlot(seurat_obj, reduction="umap", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()  
  p_cluster_tsne = DimPlot(seurat_obj, reduction="tsne", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()
  
  # save_dir = paste0(dirname(dirname(file_path)), "/seurat")
  save_dir = paste0(dirname(dirname(file_path)), "/sct")
  if(!dir.exists(save_dir)) { dir.create(save_dir, recursive = TRUE) }
  
  # output_path = file.path(save_dir, "seurat_sct.rds")
  output_path = file.path(save_dir, "seurat_sct_all_genes.rds")
  saveRDS(seurat_obj, output_path)
  
  # pdf(file.path(save_dir, "cluster.pdf"), width=7, height=6)
  # DimPlot(seurat_obj, reduction="umap", pt.size = 0.3, label = TRUE) +
  #   labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()  
  # DimPlot(seurat_obj, reduction="tsne", pt.size = 0.3, label = TRUE) +
  #   labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()
  # dev.off()

  ggsave(filename = paste0(save_dir, "/cluster_umap.png"), plot = p_cluster_umap, width=6, height=5)
  ggsave(filename = paste0(save_dir, "/cluster_tsne.png"), plot = p_cluster_tsne, width=6, height=5)

  message(paste("Processed and saved: ", save_dir))
}

# 파일 경로 탐색
get_seurat_sct = function(root_path) {
  files = dir_ls(root_path, recurse = TRUE, regexp = "seurat_sct_all_genes.rds")
  return(files)
}

# Doublet finder
doublet_filtering = function(seurat_path) {
  
  message(paste("Start with sample: ", basename(dirname(dirname(seurat_path)))))
  seurat_obj = readRDS(seurat_path)
  PCpercentage = seurat_obj@reductions$pca@stdev / sum(seurat_obj@reductions$pca@stdev) * 100
  PCcumulative = cumsum(PCpercentage)
  pc1 = which(PCcumulative > 90 & PCpercentage < 5)[1]
  pc2 = sort(which((PCpercentage[1:length(PCpercentage)-1] - PCpercentage[2:length(PCpercentage)]) > 0.1),  decreasing = T)[1] + 1
  nPC = min(pc1, pc2)
  nPC = max(nPC, 10)
  cat('selected PC value :', nPC, "\n")
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:nPC, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_pK <- bcmvn[which.max(bcmvn$BCmetric)[1], ]$pK
  nPK <- as.numeric(levels( max_pK ))[ max_pK ]
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  cells = ncol(seurat_obj)
  
  if(cells < 500 ){
    Doublet_percentage = 0.004
  } else if(500 <= cells & cells < 1000){
    Doublet_percentage = 0.008
  } else if(1000 <= cells & cells < 2000){
    Doublet_percentage = 0.016
  } else if(2000 <= cells & cells < 3000){
    Doublet_percentage = 0.024
  } else if(3000 <= cells & cells < 4000){
    Doublet_percentage = 0.032
  } else if(4000 <= cells & cells < 5000){
    Doublet_percentage = 0.040
  } else if(5000 <= cells & cells < 6000){
    Doublet_percentage = 0.048
  } else if(6000 <= cells & cells < 7000){
    Doublet_percentage = 0.056
  } else if(7000 <= cells & cells < 8000){
    Doublet_percentage = 0.064
  } else if(8000 <= cells & cells < 9000){
    Doublet_percentage = 0.072
  } else {
    Doublet_percentage = 0.080
  }

  # doublet rate : https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled

  homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)        ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(Doublet_percentage*nrow(seurat_obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:nPC, pN = 0.25, pK = nPK, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  cat('Total', cells, 'cells, Expected doublet percentage :', Doublet_percentage)
  
  n <-grepl('^pANN', colnames(seurat_obj@meta.data))
  stopifnot(sum(n)==1)
  nPANN = colnames(seurat_obj@meta.data)[n]
  nClass = paste0("DF.classifications_", substr(nPANN, 6, nchar(nPANN)))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:nPC, pN = 0.25, pK = nPK, nExp = nExp_poi.adj, reuse.pANN = nPANN, sct = TRUE)
  
  p <-grepl(nClass, colnames(seurat_obj@meta.data))
  stopifnot(sum(p)==1)
  colnames(seurat_obj@meta.data)[p] <-'DFclass'
  
  # pdf(sprintf("%s-seurat-DimPlot-doubletFinder-%s.pdf", output_barcode, as.character(nPC)))
  # DimPlot(seurat_obj, group.by="DFclass", cols=c("red", "lightgray"), raster = TRUE) + NoAxes()
  # dev.off()
  n_singlet = length(seurat_obj@meta.data$DFclass[seurat_obj@meta.data$DFclass=="Singlet"])
  n_doublet = length(seurat_obj@meta.data$DFclass[seurat_obj@meta.data$DFclass=="Doublet"])
  
  save_dir = paste0(dirname(dirname(seurat_path)), "/sct")
  if(!dir.exists(save_dir)) { dir.create(save_dir, recursive = TRUE) }
  
  sample = basename(dirname(dirname(seurat_path)))
  p_cluster_umap = DimPlot(seurat_obj, reduction="umap", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()  
  
  p_doublet = DimPlot(seurat_obj, group.by="DFclass", cols=c("red", "lightgray"),  pt.size = 0.3) + NoAxes() +
    labs(title = sample, 
         subtitle = paste0("Doublet proportion: ", Doublet_percentage, "\n",
                           "Doublet: ", n_doublet, "\n",
                           "Singlet: ", n_singlet ),
         caption = paste0("selected PC: ", nPC, "\n", 
                          "selected pK: ", nPK) )
  ggsave(filename = paste0(save_dir, "/doublet.png"), plot = p_doublet, width=10, height=5)
  
  seurat_singlet = subset(seurat_obj, subset=DFclass=="Singlet")
  output_path = file.path(save_dir, "seurat_sct_singlet.rds")
  saveRDS(seurat_singlet, output_path)
  
  sample = basename(dirname(dirname(seurat_path)))
  p_cluster_umap = DimPlot(seurat_singlet, reduction="umap", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()  
  p_cluster_tsne = DimPlot(seurat_singlet, reduction="tsne", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("selected PC: ", nPC)) + theme_classic()
  
  ggsave(filename = paste0(save_dir, "/cluster_umap_singlet.png"), plot = p_cluster_umap, width=6, height=5)
  ggsave(filename = paste0(save_dir, "/cluster_tsne_singlet.png"), plot = p_cluster_tsne, width=6, height=5)
  
  message(paste("Processed and saved: ", save_dir))
}

# 메인 함수 1
main_process = function(root_path) {
  files = get_files(root_path)
  for (file in files) {
    make_seurat_obj(file)
  }
}

# 메인 함수 2
doublet_finder_process = function(root_path) {
  files = get_seurat_sct(root_path)
  for (file in files) {
    doublet_filtering(file)
  }
}

# 실행
root_path = "/Volumes/WindySSD/YUGBM_OC/cellranger/run_cellranger_count/"
main_process(root_path)
doublet_finder_process(root_path)

# 메인 함수 2 실행
#root_path = "/Volumes/WindySSD/YUGBM_OC/cellranger/run_cellranger_count_newsample/CN11-C"






# 합친거









# 





# read_file = function(f){
#   dat.fl = readRDS(f)}
# 
# list = list.files(path="/Volumes/WindySSD/YUGBM_OC/cellranger/run_cellranger_count/", recursive=T, pattern="seurat_sct_singlet.rds", full.names=T)
# 
# data = sapply(list, read_file)
# 
# sample_code = c("CN10-V", "CN9-V", 
#                 "GC14-C", "GC14-T", 
#                 "GN1-V", "GN1-T", 
#                 "GN15-T1", "GN15-T2", "GN15-V", 
#                 "GN16-T1", "GN16-T2", "GN16-V",
#                 "GN17-T1", "GN17-T2", "GN17-V", 
#                 "GN2-V", "GN2-T", 
#                 "GN3-T", "GN3-V", 
#                 "GN5-V", "GN5-T", 
#                 "GN6-V", "GN6-T", 
#                 "GN7-C", "GN7-V", "GN7-T", 
#                 "GN8-T", "GN8-V")
# 
# data_merged = merge(x = data[[1]],
#                     y = data[2:length(data)],
#                     add.cell.ids = sample_code)
# saveRDS(data_merged, "/Users/home/Desktop/project/gbmoc/data/merge/sct_merged_v1.RDS")
# 
# 

# # Version 2 test (ver2.1에서 GN6 sample 제외)
# list = list.files(path="/Volumes/WindySSD/YUGBM_OC/cellranger/ver2", recursive=T, pattern="seurat_sct_singlet.rds", full.names=T)
# read_file = function(f){
#   dat.fl = readRDS(f)}
# data = sapply(list, read_file)
# 
# sample_code = c("CN10-V", "CN9-V", 
#                 "GC14-C", "GC14-T", 
#                 "GN1-V", "GN1-T", 
#                 "GN15-T1", "GN15-V", 
#                 "GN16-T1", "GN16-V",
#                 "GN17-T1", "GN17-V", 
#                 "GN2-V", "GN2-T", 
#                 "GN3-T", "GN3-V", 
#                 "GN7-C", "GN7-V", "GN7-T", 
#                 "GN8-T", "GN8-V")
# 
# merged = merge(x = data[[1]],
#                     y = data[2:length(data)],
#                     add.cell.ids = sample_code)
# saveRDS(merged, "/Users/home/Desktop/project/gbmoc/data/merge/sct_merged_v2.RDS")
# merged = readRDS("/Users/home/Desktop/project/gbmoc/data/merge/sct_merged_v2.RDS")
# 
# # merged[["RNA"]] = split(merged[["RNA"]], f=merged$sample_barcode)
# merged = readRDS("/Users/home/Desktop/project/gbmoc/data/merge/sct_merged_v2.RDS")
# merged = SCTransform(merged, vst.flavor="v2")
# merged = RunPCA(merged)
# saveRDS(merged, "/Users/home/Desktop/project/gbmoc/data/merge/sct_merged_pca_v2.RDS")
# 
# # Integrate (Harmony)
# integrated = IntegrateLayers(
#   object = merged, method = HarmonyIntegration,
#   orig.reduction = "pca", new.reduction = "harmony",
#   normalization.method = "SCT", 
#   assay = "SCT", verbose = TRUE)

# Version 2.3 : 23 samples
# Version 2.4 : 24 samples
pkgs = c("Seurat", "sctransform", "glmGamPoi",  
         "dplyr", "plyr", "tibble", "fs", "ggplot2", "parallel", "future")
lapply(pkgs, library, character.only = TRUE)
options(future.globals.maxSize = 50*1024^3)

read_file = function(f){
  dat.fl = readRDS(f)}
version_num = "v2.4"
list = list.files(path=sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s",version_num), recursive=T, pattern="seurat_sct_singlet.rds", full.names=T)
data = sapply(list, read_file)

sample_code = c("CN10-V", "CN10-C", "CN11-C", "CN9-V", 
                "GC14-C",
                "GN1-V", "GN1-T", 
                "GN15-T1", "GN15-V", 
                "GN16-T1", "GN16-V",
                "GN17-T1", "GN17-V", 
                "GN2-V", "GN2-T", 
                "GN3-T", "GN3-V", 
                "GN6-V", "GN6-T",
                "GN7-C", "GN7-V", "GN7-T", 
                "GN8-T", "GN8-V")

merged = merge(x = data[[1]],
               y = data[2:length(data)],
               add.cell.ids = sample_code)

merged@meta.data$sample_barcode = sapply(strsplit(rownames(merged@meta.data), "_"), `[`, 1)
merged@meta.data$patient_barcode = sapply(strsplit(merged@meta.data$sample_barcode, "-"), `[`, 1)
print( sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_merged_%s.RDS", version_num, version_num) )
saveRDS(merged, sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_merged_%s.RDS", version_num, version_num))

merged = readRDS(sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_merged_%s.RDS", version_num, version_num))
merged = SCTransform(merged, vst.flavor="v2")
merged = RunPCA(merged)
saveRDS(merged, sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_merged_%s.RDS", version_num, version_num))


integrated = IntegrateLayers(
  object = merged, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  normalization.method = "SCT", 
  assay = "SCT", verbose = TRUE)
# saveRDS(integrated, sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_integrated_%s.RDS", version_num, version_num))

# Select PC value
ElbowPlot(integrated)
PCpercentage = integrated@reductions$pca@stdev / sum(integrated@reductions$pca@stdev) * 100
PCcumulative = cumsum(PCpercentage)
pc1 = which(PCcumulative > 90 & PCpercentage < 5)[1]
pc2 = sort(which((PCpercentage[1:length(PCpercentage)-1] - PCpercentage[2:length(PCpercentage)]) > 0.1),  decreasing = T)[1] + 1
nPC = min(pc1, pc2)
cat('selected PC value :', nPC)

integrated = FindNeighbors(integrated, reduction = "harmony", dims=1:nPC)
integrated = FindClusters(integrated, resolution = 0.5, cluster.name = "harmony_clusters")
integrated = RunUMAP(integrated, reduction = "harmony", dims=1:nPC, resolution=0.5, reduction.name = "umap.harmony")
integrated = RunTSNE(integrated, reduction = "harmony", dims=1:nPC, resolution=0.5, reduction.name = "tsne.harmony")

integrated@meta.data$sample_barcode = sapply(strsplit(rownames(integrated@meta.data), "_"), `[`, 1)
integrated@meta.data$patient_barcode = sapply(strsplit(integrated@meta.data$sample_barcode, "-"), `[`, 1)
integrated@meta.data$sample_type1 = sapply(strsplit(integrated@meta.data$sample_barcode, "-"), `[`, 2)

DimPlot(integrated, reduction="tsne.harmony")
DimPlot(integrated, reduction="umap.harmony")
DimPlot(integrated, reduction="tsne.harmony", group.by="sample_barcode")
DimPlot(integrated, reduction="umap.harmony", group.by="sample_barcode")
DimPlot(integrated, reduction="tsne.harmony", group.by="patient_barcode")
DimPlot(integrated, reduction="umap.harmony", group.by="patient_barcode")
cat(sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_integrated_%s.RDS", version_num, version_num))
saveRDS(integrated, sprintf("/Volumes/WindySSD/YUGBM_OC/cellranger/%s/sct_integrated_%s.RDS", version_num, version_num))




