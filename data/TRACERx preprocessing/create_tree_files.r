#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("openxlsx", dependencies = TRUE)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: ./create_tree_files.r <input_path> <output_path>")
}

input_path <- args[1]
outdir     <- args[2]

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Please replace /path/to with the path to the TRACERx data
data <- readRDS(file.path(input_path, "data/20221109_TRACERx421_phylogenetic_trees.rds"))

mutation_data <- readRDS(file.path(input_path, "data/LCCE_converted_data/20221109_TRACERx421_mutation_table.rds"))
#print(mutation_data[mutation_data$patient_id == "CRUK0012", ])

patient_data <- readRDS(file.path(input_path, "data/20221109_TRACERx421_all_patient_df.rds"))
# Corrected tree is the default tree reported by CONIPHER
#print(data[["CRUK0386"]]$graph_pyclone$Corrected_tree)
# In alt_trees are all trees reported by CONIPHER
#print(data[["CRUK0386"]]$graph_pyclone$alt_trees)

keys <- names(data)
for (i in seq_len(length(keys))) {
  print(i)
  key <- keys[i]
  # Uncomment the two lines below if you want to extract the default trees and comment out the two lines after that
  #file_name <- file.path(outdir, paste0(key, '.xlsx'))
  #write.xlsx(data[[key]]$graph_pyclone$Corrected_tree, file_name)
  file_name <- file.path(outdir, paste0(key, '.xlsx'))
  write.xlsx(data[[key]]$graph_pyclone$alt_trees, file_name)
}