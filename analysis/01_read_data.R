
## 01: Load data from cellranger output into Seurat object and save on disk

# load libs
library(Seurat)
library(scUtils)
library(dplyr)

data_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq/28022023/"
feature_column_idx = 2
project_name = "AgrpPomcNTS"

# load sample names
# sampleNames_tab = data.table::fread(paste0(data_path,"Sample_Names.tab"),data.table = FALSE)
# sample_table = sampleNames_tab %>% dplyr::rename(SID = V1, SampleName = `CCG Sample ID    Sample Name`) %>% dplyr::mutate(SID = paste0("SID",SID))

# load metadata
sample_metadata = data.table::fread(paste0(data_path,"2023_03_almu_sample_metadata.txt"),data.table = F)
colnames(sample_metadata)[1] = "SID"
sample_ids = paste0("SID",sample_metadata$SID)[sample_metadata$Treatment %in% c("WT_CNO","AgRP-Gq_POMC-Gi_CNO")]

# load files
all_files = list.files(path = data_path,recursive = TRUE,full.names = TRUE)
all_sample_seurats=list()
for(sample_id in sample_ids){
  message(sample_id)
  sample_files = all_files[grepl(sample_id,all_files)]
  sample_files = sample_files[grepl("filtered_feature_bc_matrix",sample_files)]
  sample_counts <- scUtils::Read10xFormat(mtx = sample_files[grepl("matrix.mtx",sample_files)],
                                          cells = sample_files[grepl("barcodes",sample_files)],
                                          features = sample_files[grepl("features",sample_files)],
                                          feature.column = feature_column_idx)
  # make seurat
  # add run name to column names
  colnames(sample_counts) = paste0(colnames(sample_counts),"_",sample_id)
  
  # make seurat object
  sample_seurat = SeuratObject::CreateSeuratObject(sample_counts,project = sample_id,min.cells = 0,  min.features = 0 )
  sample_seurat@meta.data$Cell_ID = colnames(sample_seurat)
  sample_seurat@meta.data$Sample_ID = sample_id
  sample_seurat@meta.data$MouseLine = sample_metadata$Line[sample_metadata$SID == gsub("SID","",sample_id)]
  sample_seurat@meta.data$Condition = sample_metadata$Treatment[sample_metadata$SID == gsub("SID","",sample_id)]
  all_sample_seurats[[sample_id]] = sample_seurat
}
# merge seurat objects:
dataset_seurat <- merge(all_sample_seurats[[1]], y = all_sample_seurats[2:length(all_sample_seurats)], project = project_name)
# add mt
dataset_seurat[["percent.mt"]] <- Seurat::PercentageFeatureSet(dataset_seurat, pattern = "^mt-")
dataset_seurat@meta.data$Condition2 = dataset_seurat@meta.data$Condition
dataset_seurat@meta.data$Condition2[dataset_seurat@meta.data$Condition2 == "AgRP-Gq_POMC-Gi_CNO"] = "AgrpGq_POMCGi"

## save
outputdata_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/"
dir.create(outputdata_path)
saveRDS(dataset_seurat, paste0(outputdata_path,"almu_2023_04_AgrpPomc_raw.rds"))


