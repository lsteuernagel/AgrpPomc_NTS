## 06: Prepare files and objects for GEO supplementary data upload and publication source data

##########
### Load data
##########

## libs
library(tidyverse)
library(Seurat)

## project settings with local file paths
raw_file_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/","almu_2023_04_AgrpPomc_raw.rds")
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"
figure_path = "paper_figures/"
source_path = "paper_sourcedata/"

# seurat
AgrpPomcNTS_withIEG = readRDS(paste0(project_path,"AgrpPomcNTS_withIEG.rds"))
# degs
degs_all = data.table::fread(paste0(project_path,"deg_nebula_all.txt"),data.table = F)
degs_all = degs_all %>% dplyr::arrange(cluster,pval_adj,pval)
degs_all$cluster = as.character(degs_all$cluster)


##########
### Export supplementary data for GEO
##########

## saving large files for GEO upload outside git repo!


##########
### Export source data as xlsx
##########

all_source_files = list.files(source_path,pattern = "source_",full.names = TRUE) 
all_source_files_short = list.files(source_path,pattern = "source_",full.names = FALSE) 
source_table_list = list()
source_table_list_ext = list()
for(i in 1:length(all_source_files)){
  current_table = data.table::fread(paste0(all_source_files[i]),data.table = FALSE) 
  current_name = gsub("source_|\\.txt","",all_source_files_short[i])
  if(grepl("ext",all_source_files_short[i])){
    source_table_list_ext[[current_name]] = current_table
  }else{
    source_table_list[[current_name]] = current_table
  }
}

# save as xlsx
WriteXLS::WriteXLS(x = source_table_list,ExcelFileName = paste0(source_path,"Figure_5_source.xlsx"),col.names=TRUE)
WriteXLS::WriteXLS(x = source_table_list_ext,ExcelFileName = paste0(source_path,"ExtFigure_10_source.xlsx"),col.names=TRUE)

##########
### Also save R packages
##########

## load all relevant functions and scripts
session_info = sessionInfo() 
session_info$R.version$version.string
writeLines(capture.output(sessionInfo()), "paper_sourcedata/sessionInfo.txt")
# I just copy the main loaded packages (not all dependcies via namespace) and cat them in a second file for the reporting summary etc.:
all_packages = unlist(strsplit(
  gsub("\\[[0-9]*\\]","","[1] ggh4x_0.2.4                 nebula_1.2.1                scDblFinder_1.12.0          SingleCellExperiment_1.20.0 SummarizedExperiment_1.28.0
 [6] Biobase_2.58.0              GenomicRanges_1.50.2        GenomeInfoDb_1.34.4         IRanges_2.32.0              S4Vectors_0.36.1           
[11] BiocGenerics_0.44.0         MatrixGenerics_1.10.0       matrixStats_0.63.0-9003     Matrix_1.5-3                scUtils_0.0.1              
[16] magrittr_2.0.3              SeuratObject_4.1.3          Seurat_4.3.0                forcats_0.5.2               stringr_1.5.0              
[21] dplyr_1.1.0                 purrr_1.0.1                 readr_2.1.3                 tidyr_1.3.0                 tibble_3.1.8               
[26] tidyverse_1.3.2             ggpubr_0.6.0                ggplot2_3.4.2         "
  ),split=" "
))
all_packages = all_packages[all_packages!="" & all_packages != "\n"]
all_packages_df = data.frame(package = gsub("_"," - version:",all_packages))
data.table::fwrite(all_packages_df,file = "paper_sourcedata/R_packages.txt",sep="\t")



