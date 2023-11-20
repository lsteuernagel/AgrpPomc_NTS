## 06: Preepare files and objects for GEO supplementalry data upload

## libs
library(tidyverse)
library(Seurat)

## project settings
raw_file_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/","almu_2023_04_AgrpPomc_raw.rds")
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"

# seurat
AgrpPomcNTS_withIEG = readRDS(paste0(project_path,"AgrpPomcNTS_withIEG.rds"))
# degs
degs_all = data.table::fread(paste0(project_path,"deg_nebula_all.txt"),data.table = F)
degs_all = degs_all %>% dplyr::arrange(cluster,pval_adj,pval)
degs_all$cluster = as.character(degs_all$cluster)
