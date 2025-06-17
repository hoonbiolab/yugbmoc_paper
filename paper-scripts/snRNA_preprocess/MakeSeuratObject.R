###########################################################################################
# Packages
###########################################################################################
pkgs = c("Seurat", "dplyr", "plyr", "tibble", "fs", "ggplot2", "DoubletFinder")
lapply(pkgs, library, character.only = TRUE)

###########################################################################################
# Functions
###########################################################################################
read_file = function(f){
  dat.fl = readRDS(f)}

get_seurat_file_path = function(root_path) {
  file_path = dir_ls(root_path, recurse = TRUE, regexp = "seurat_singlet.rds$")
  return(file_path)}

get_h5_files = function(root_path) {
  files = dir_ls(root_path, recurse = TRUE, regexp = "filtered_feature_bc_matrix.h5$")
  return(files)}

get_infercnv_obj_path = function(infercnv_outpath) {
  file_path = dir_ls(infercnv_outpath, recurse = TRUE, regexp = "run.final.infercnv_obj$")
  return(dirname(file_path))}

###########################################################################################
# QC threshold
###########################################################################################
min_cells = 10
min_features = 200
max_percent_mt = 5
min_nFeature_RNA = 200
max_nFeature_RNA = 10000
min_nCount_RNA = 1000
max_nCount_RNA = 15000

###########################################################################################
# 10x -> Seurat object : Preprocess (EACH SAMPLE)
###########################################################################################
preprocessed_seurat_obj = function(file_path) {
  message(paste("Start with sample: ", basename(dirname(dirname(file_path))), "\n"))
  seurat_obj = CreateSeuratObject(counts = Read10X_h5(file_path, use.names = TRUE),
                                  min.cells = min_cells, min.features = min_features,
                                  project = "YUGBM")
  seurat_obj[["percent.mt"]] = PercentageFeatureSet(seurat_obj, pattern="^MT-")
  
  p_before_filtering = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0) +
    labs(caption = "Before filtering")
  
  seurat_obj <- subset(seurat_obj, 
                       subset = 
                         percent.mt < max_percent_mt &
                         nFeature_RNA > min_nFeature_RNA & 
                         nFeature_RNA < max_nFeature_RNA &
                         nCount_RNA > min_nCount_RNA & 
                         nCount_RNA < max_nCount_RNA )
  
  p_after_filtering = VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, alpha = 0) +
    labs(caption = "After filtering")
  
  p_filtering = p_before_filtering / p_after_filtering
  
  seurat_obj = NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  seurat_obj = FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
  seurat_obj = ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj = RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
  # Normalized values = seurat_obj[["RNA"]]$data
  # Scale data = seurat_obj[["RNA"]]$scale.data : mean 0 variance 1
  
  PCpercentage = seurat_obj@reductions$pca@stdev / sum(seurat_obj@reductions$pca@stdev) * 100
  PCcumulative = cumsum(PCpercentage)
  pc1 = which(PCcumulative > 90 & PCpercentage < 5)[1]
  pc2 = sort(which((PCpercentage[1:length(PCpercentage)-1] - PCpercentage[2:length(PCpercentage)]) > 0.05),  decreasing = T)[1] + 1
  nPC = min(pc1, pc2)
  cat('selected PC value :', nPC, "\n")

  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list <- paramSweep(seurat_obj, PCs = 1:nPC, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_pK <- bcmvn[which.max(bcmvn$BCmetric)[1], ]$pK
  nPK <- as.numeric(levels( max_pK ))[ max_pK ]
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  cells = ncol(seurat_obj)
  
  if(cells < 500 ){
    Doublet_proportion = 0.004
  } else if(500 <= cells & cells < 1000){
    Doublet_proportion = 0.008
  } else if(1000 <= cells & cells < 2000){
    Doublet_proportion = 0.016
  } else if(2000 <= cells & cells < 3000){
    Doublet_proportion = 0.024
  } else if(3000 <= cells & cells < 4000){
    Doublet_proportion = 0.032
  } else if(4000 <= cells & cells < 5000){
    Doublet_proportion = 0.040
  } else if(5000 <= cells & cells < 6000){
    Doublet_proportion = 0.048
  } else if(6000 <= cells & cells < 7000){
    Doublet_proportion = 0.056
  } else if(7000 <= cells & cells < 8000){
    Doublet_proportion = 0.064
  } else if(8000 <= cells & cells < 9000){
    Doublet_proportion = 0.072
  } else {
    Doublet_proportion = 0.080
  }
  # doublet rate : https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled
  
  homotypic.prop <- modelHomotypic(seurat_obj@meta.data$seurat_clusters)        ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(Doublet_proportion*nrow(seurat_obj@meta.data))              ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:nPC, pN = 0.25, pK = nPK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  cat('Total', cells, 'cells, Expected doublet proportion :', Doublet_proportion)
  
  n <-grepl('^pANN', colnames(seurat_obj@meta.data))
  stopifnot(sum(n)==1)
  nPANN = colnames(seurat_obj@meta.data)[n]
  nClass = paste0("DF.classifications_", substr(nPANN, 6, nchar(nPANN)))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seurat_obj <- doubletFinder(seurat_obj, PCs = 1:nPC, pN = 0.25, pK = nPK, nExp = nExp_poi.adj, reuse.pANN = nPANN, sct = FALSE)
  
  p <-grepl(nClass, colnames(seurat_obj@meta.data))
  stopifnot(sum(p)==1)
  colnames(seurat_obj@meta.data)[p] <-'DFclass'
  
  n_singlet = length(seurat_obj@meta.data$DFclass[seurat_obj@meta.data$DFclass=="Singlet"])
  n_doublet = length(seurat_obj@meta.data$DFclass[seurat_obj@meta.data$DFclass=="Doublet"])
  
  save_dir = paste0(dirname(dirname(file_path)), "/RDS")
  if(!dir.exists(save_dir)) { dir.create(save_dir, recursive = TRUE) }
  plot_dir = paste0(dirname(dirname(file_path)), "/Plot")
  if(!dir.exists(plot_dir)) { dir.create(plot_dir, recursive = TRUE) }
  
  sample = basename(dirname(dirname(file_path)))
  
  seurat_obj = RunUMAP(seurat_obj, dims = 1:nPC, verbose = FALSE)
  seurat_obj = RunTSNE(seurat_obj, dims = 1:nPC, verbose = FALSE)
  seurat_obj = FindNeighbors(seurat_obj, dims = 1:nPC, verbose = FALSE)
  seurat_obj = FindClusters(seurat_obj, verbose = FALSE)
  
  p_doublet = DimPlot(seurat_obj, group.by="DFclass", cols=c("red", "lightgray"),  pt.size = 0.3) + NoAxes() +
    labs(title = sample, 
         subtitle = paste0("Doublet proportion: ", Doublet_proportion, "\n",
                           "Singlet: ", n_singlet, "\n",
                           "Doublet: ", n_doublet),
         caption = paste0("selected PC: ", nPC, "\n", 
                          "selected pK: ", nPK) )
  
  seurat_obj = subset(seurat_obj, subset=DFclass=="Singlet")
  seurat_obj@meta.data <- seurat_obj@meta.data[, !grepl("^DF\\.classifications", colnames(seurat_obj@meta.data))]
  seurat_obj@meta.data <- seurat_obj@meta.data[, !grepl("^pANN", colnames(seurat_obj@meta.data))]
  
  p_cluster_umap = DimPlot(seurat_obj, reduction="umap", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("cells: ", ncol(seurat_obj), "\n",
                                           "selected PC: ", nPC)) + theme_classic()  
  p_cluster_tsne = DimPlot(seurat_obj, reduction="tsne", pt.size = 0.3, label = TRUE) +
    labs(title = sample, subtitle = paste0("cells: ", ncol(seurat_obj), "\n",
                                           "selected PC: ", nPC)) + theme_classic()  
  
  ggsave(filename = paste0(plot_dir, "/qc.png"), plot = p_filtering, width=6, height = 8)
  ggsave(filename = paste0(plot_dir, "/doublet.png"), plot = p_doublet, width=10, height=5)
  ggsave(filename = paste0(plot_dir, "/cluster_umap.png"), plot = p_cluster_umap, width=6, height=5)
  ggsave(filename = paste0(plot_dir, "/cluster_tsne.png"), plot = p_cluster_tsne, width=6, height=5)
  
  seurat_obj_resize = CreateSeuratObject(counts=as.matrix(seurat_obj[["RNA"]]$counts), project=)
  output_path = file.path(save_dir, "seurat_singlet.rds")
  saveRDS(seurat_obj_resize, output_path)
  message(paste("Processed Finish! Saved in: ", save_dir))
  rm(seurat_obj)
  rm(seurat_obj_resize)
}

preprocess = function(root_path) {
  files = get_h5_files(root_path)
  for (file in files) {
    preprocessed_seurat_obj(file)
  }
}

