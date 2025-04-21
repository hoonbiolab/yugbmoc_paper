# Load necessary libraries
library(dplyr)
library(plyr)
library(tidyverse)
library(RNAseqCNV)

# Step 1: Process RSEM gene expression files
# Get the list of RSEM gene result files
rsemfile_list <- list.files('~/Dropbox/HoonLab/YUMED_RNA/workflow_out/RSEM_v1.3.3', 
                            pattern = '.rsem.genes.results.gz$', 
                            recursive = TRUE, 
                            full.names = TRUE)

# Load and filter RSEM files
Rsem_df_filtered <- ldply(rsemfile_list, function(file_path) {
  cat("Processing:", file_path, "\n")
  d <- read.table(gzfile(file_path), sep = "\t", header = TRUE)
  # Select gene_id and expected_count columns
  d %>% select(gene_id, expected_count) %>% 
    mutate(fname = basename(file_path))
})

# Add sample name by removing file extension
Rsem_df_filtered <- Rsem_df_filtered %>% 
  mutate(Sample_Name = gsub("\\.rsem.genes.results.gz", "", fname))

# Convert expected_count to integer
Rsem_df_filtered$expected_count <- as.integer(Rsem_df_filtered$expected_count)

# Display the first row to verify
head(Rsem_df_filtered, 1)

# Step 2: Write count files for each sample
# Set working directory
setwd("~/Dropbox/HoonLab/YUMED_RNA/workflow_out/RNAseqCNV/")

# Create count files for each sample
sample_name_list <- unique(Rsem_df_filtered$Sample_Name)
for (sample_name in sample_name_list) {
  sample_data <- Rsem_df_filtered %>% 
    filter(Sample_Name == sample_name) %>% 
    select(gene_id, expected_count)
  
  # Write each sample's data to a text file
  write.table(sample_data, 
              file = paste0("count_files/", sample_name, ".txt"), 
              quote = FALSE, sep = '\t', 
              row.names = FALSE, col.names = FALSE)
}

# Step 3: Process VCF files
# Get the list of VCF files
vcffile_list <- list.files('~/Dropbox/HoonLab/YUMED_RNA/workflow_out/gatk_germline_indel/variant_filtered_vcfs/', 
                           pattern = '.vcf.gz$', 
                           recursive = TRUE, 
                           full.names = TRUE)

# Extract sample names from VCF file names
vcf_sample_list <- gsub(".variant_filtered.vcf.gz", "", basename(vcffile_list))

# Step 4: Generate count files for samples that have VCF data
for (sample_name in vcf_sample_list) {
  sample_data <- Rsem_df_filtered %>% 
    filter(Sample_Name == sample_name) %>% 
    select(gene_id, expected_count)
  
  # Write each sample's data to a text file
  write.table(sample_data, 
              file = paste0("count_files/", sample_name, ".txt"), 
              quote = FALSE, sep = '\t', 
              row.names = FALSE, col.names = FALSE)
}

# Filter the RSEM data to include only samples with VCF files
Rsem_df_filtered <- Rsem_df_filtered %>% 
  filter(Sample_Name %in% vcf_sample_list)

length(unique(Rsem_df_filtered$Sample_Name))


# Step 5: Create metadata for RNAseqCNV
# Match the count files and VCF files by sample name
meta_data_df <- data.frame(Sample_Name = vcf_sample_list, 
                           Count_File = paste0(vcf_sample_list, ".txt"), 
                           Vcf_File = basename(vcffile_list))

# Write the metadata to a file
write.table(meta_data_df, 
            file = "metadata/YUMED_meta.txt", 
            quote = FALSE, sep = '\t', 
            row.names = FALSE, col.names = FALSE)

# Step 6: Run RNAseqCNV
# Execute RNAseqCNV with the generated config and metadata files
RNAseqCNV_wrapper(config = "config/RNAseqCNV.config", 
                  metadata = "YUMED_meta.txt", 
                  snv_format = "vcf", 
                  genome_version = "hg19", 
                  CNV_matrix = TRUE)
