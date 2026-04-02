# Supplementary Fig. 6C

# Load necessary libraries
pkgs = c('dplyr', 'readr', 'tidyr', 'reshape2', 'data.table', 'tidyverse', 'tibble',
         'ggplot2', 'pheatmap', 'ComplexHeatmap', 'scales', 'colorRamp2')
lapply( pkgs, library, character.only=T )

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Step 1
samples = readxl::read_excel("~/Projects/gbm-oc/RNA/RNAseqCNV/manifest/gbm-oc-rna-datasets.xlsx") %>%
  dplyr::select(c("aliquot_barcode", "patient_barcode", "patient_barcode_legacy", "sample_type", "stage", "gender_legacy"))

rnacnv = read.table("~/Projects/gbm-oc/RNA/RNAseqCNV/data/log2_fold_change_per_arm.tsv", sep='\t', header=T) %>%
  mutate(chr = paste0(chr, arm)) %>% 
  select(-c("arm")) %>%
  remove_rownames %>%
  column_to_rownames(var="chr")

rnacnv.t = as.data.frame(t(rnacnv[, -1])) %>%
  mutate(aliquot_barcode = rownames(.)) %>%
  mutate(aliquot_barcode = gsub("[.]", "-", aliquot_barcode))

samples.rnacnv = left_join(samples, rnacnv.t, join_by(aliquot_barcode == aliquot_barcode)) %>%
  na.omit() %>%
  select(-c("patient_barcode_legacy", "stage", "gender_legacy", "Xp", "Xq")) %>%
  mutate(stage = substr(aliquot_barcode, 20, 20)) %>%
  filter(stage != "R") 


# Step 2

# SVZ-control
rnacnv.mat.svzc = samples.rnacnv %>% column_to_rownames(var="aliquot_barcode") %>%
  filter(sample_type =="SVZ-control") %>%
  select(-c('sample_type', "patient_barcode", "stage")) %>%
  data.matrix(., rownames.force=T) %>%
  t(.)

rnacnv.annot.svzc = samples.rnacnv %>%
  filter(sample_type =="SVZ-control") %>%
  select(c("aliquot_barcode", "sample_type", "patient_barcode"))

rnacnv.annot.svzc = HeatmapAnnotation(
  sample = rnacnv.annot.svzc$sample_type,
  col = list(sample = c("SVZ" = "orange", "GBM" = "red", "SVZ-control" = "white")))

# GBM
rnacnv.mat.gbm = samples.rnacnv %>% column_to_rownames(var="aliquot_barcode") %>%
  filter(sample_type =="GBM") %>%
  select(-c('sample_type', "patient_barcode", "stage")) %>%
  data.matrix(., rownames.force=T) %>%
  t(.)

rnacnv.annot.gbm = samples.rnacnv %>%
  filter(sample_type =="GBM") %>%
  select(c("aliquot_barcode", "sample_type", "patient_barcode"))

rnacnv.annot.gbm = HeatmapAnnotation(
  sample = rnacnv.annot.gbm$sample_type,
  col = list(sample = c("SVZ" = "orange", "GBM" = "red", "SVZ-control" = "white")))

# SVZ
rnacnv.mat.svz = samples.rnacnv %>% column_to_rownames(var="aliquot_barcode") %>%
  filter(sample_type =="SVZ") %>%
  select(-c('sample_type', "patient_barcode", "stage")) %>%
  data.matrix(., rownames.force=T) %>%
  t(.)

rnacnv.annot.svz = samples.rnacnv %>%
  filter(sample_type =="SVZ") %>%
  select(c("aliquot_barcode", "sample_type", "patient_barcode"))

rnacnv.annot.svz = HeatmapAnnotation(
  sample = rnacnv.annot.svz$sample_type,
  col = list(sample = c("SVZ" = "orange", "GBM" = "red", "SVZ-control" = "white")))

# Cortex
rnacnv.mat.cor = samples.rnacnv %>% column_to_rownames(var="aliquot_barcode") %>%
  filter(sample_type =="Cortex") %>%
  select(-c('sample_type', "patient_barcode", "stage")) %>%
  data.matrix(., rownames.force=T) %>%
  t(.)

rnacnv.mat.cor.average = rowSums(rnacnv.mat.cor) / ncol(rnacnv.mat.cor)
samples.rnacnv$sample_type = factor(samples.rnacnv$sample_type, levels=c("Cortex", "SVZ-control", "SVZ", "GBM"))
samples.rnacnv = samples.rnacnv %>% arrange(sample_type, `7p`)
my_col = colorRamp2(c(-0.2, 0, 0.2), c("blue", "white", "red"))

# Baseline as cortex average
rnacnv.mat.svz.adj = rnacnv.mat.svz - rnacnv.mat.cor.average
rnacnv.mat.gbm.adj = rnacnv.mat.gbm - rnacnv.mat.cor.average
rnacnv.mat.cor.adj = rnacnv.mat.cor - rnacnv.mat.cor.average
rnacnv.mat.svzc.adj = rnacnv.mat.svzc -rnacnv.mat.cor.average

hm.gbm.adj = Heatmap(rnacnv.mat.gbm.adj, 
                     column_title = "GBM",
                     cluster_rows = FALSE,
                     # cluster_columns = FALSE,
                     clustering_distance_columns = "spearman",
                     #column_km = 2, column_km_repeats = 100,
                     show_column_names = FALSE,
                     col = my_col, 
                     row_names_side = "left",
                     top_annotation = rnacnv.annot.gbm,
                     name = "CN range")

hm.svz.adj = Heatmap(rnacnv.mat.svz.adj, 
                     column_title = "SVZ",
                     cluster_rows = FALSE,
                     # cluster_columns = FALSE,
                     clustering_distance_columns = "spearman",
                     #column_km = 2, column_km_repeats = 100,
                     show_column_names = FALSE,
                     col = my_col, 
                     row_names_side = "left",
                     top_annotation = rnacnv.annot.svz,
                     show_heatmap_legend = FALSE)

hm.svzc.adj = Heatmap(rnacnv.mat.svzc.adj, 
                      column_title = "SVZ-Control",
                      cluster_rows = FALSE,
                      # cluster_columns = FALSE,
                      clustering_distance_columns = "spearman",
                      show_column_names = FALSE,
                      col = my_col, 
                      row_names_side = "left",
                      top_annotation = rnacnv.annot.svzc,
                      show_heatmap_legend = FALSE)

# Final Figure
hm_adj_list = hm.svzc.adj  + hm.svz.adj  + hm.gbm.adj   

# Save
pdf("~/Projects/gbm-oc/RNA/RNAseqCNV/plot/Heatmap-RNAseqCNV.pdf")
draw(hm_adj_list, ht_gap = unit(1, "mm"))
dev.off()
