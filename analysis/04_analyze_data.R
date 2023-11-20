## 04: Analyze and annotate data

# load libs
library(Seurat)
library(scUtils)
library(tidyverse)

# need
path_to_pipeline = "/beegfs/scratch/bruening_scratch/lsteuernagel/projects/scseq_pipeline/" # path_to_pipeline = ""

## project settings
raw_file_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/","almu_2023_04_AgrpPomc_raw.rds")
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"

## extra colors:
short_palette = as.character(palette.colors(palette = "Okabe-Ito"))
short_palette = short_palette[!short_palette %in% c("#999999","#000000")]
getOkabeItoPalette = colorRampPalette(short_palette)

## read
AgrpPomcNTS_processed = readRDS(paste0(project_path,"AgrpPomcNTS_processed.rds"))

### load markers
marker_genes = AgrpPomcNTS_processed@misc$marker_genes %>% dplyr::filter(p_val_adj < 0.0001) %>%
  dplyr::mutate(specificity = avg_log2FC * ((pct.1+0.01) / (pct.2+0.01)))

## manually update some annotation
AgrpPomcNTS_processed@meta.data$celltype_annotation[AgrpPomcNTS_processed@meta.data$preliminary_clusters  %in% c("23","29","34")] = "Choroid plexus" #  Otx2+,Folr1,Prlr+ #https://www.nature.com/articles/s41380-021-01416-3
AgrpPomcNTS_processed@meta.data$celltype_annotation[AgrpPomcNTS_processed@meta.data$preliminary_clusters %in% c("61","141")] = "Tanycyte-like"
AgrpPomcNTS_processed@meta.data$celltype_annotation[AgrpPomcNTS_processed@meta.data$preliminary_clusters %in% c("19","11")] = "low-quality-neurons"
AgrpPomcNTS_processed@meta.data$celltype_annotation[AgrpPomcNTS_processed@meta.data$celltype_annotation == "Dcn.Fibroblasts"] = "Fibroblasts"

##########
### Initial plots
##########

DimPlot(AgrpPomcNTS_processed,group.by = "celltype_annotation",label=TRUE)+NoLegend()
DimPlot(AgrpPomcNTS_processed,group.by = "SNN_scvi_native_res.20",label=TRUE,label.size = 2)+NoLegend()
DimPlot(AgrpPomcNTS_processed,group.by = "preliminary_clusters",label=TRUE,label.size = 2)+NoLegend()
DimPlot(AgrpPomcNTS_processed,group.by = "Sample_ID",shuffle = TRUE)

FeaturePlot(AgrpPomcNTS_processed,"Prlr",order=TRUE)
table(AgrpPomcNTS_processed@meta.data$Condition,AgrpPomcNTS_processed@meta.data$Sample_ID)


cellsh = AgrpPomcNTS_processed@meta.data$[AgrpPomcNTS_processed@meta.data$preliminary_clusters == "17"]
DimPlot(AgrpPomcNTS_processed,cells.highlight = cellsh)+NoLegend()

DimPlot(AgrpPomcNTS_processed,group.by = "celltype_annotation",label = TRUE)+NoLegend()

##########
### Add cluster anno
##########

## top 2 marker genes
top_2_markers = marker_genes %>%
  dplyr::group_by(cluster) %>%
  dplyr::filter(specificity > 2 & p_val_adj < 0.01) %>%
  dplyr::filter(!grepl("Rik|Gm|D7",gene)) %>%
  dplyr::slice_max(order_by = specificity,n = 2,with_ties = F) %>%
  dplyr::summarise(top_markers = paste0(gene,collapse = "/"))

top_2_markers$cluster_annotation = paste0(top_2_markers$cluster,": ",top_2_markers$top_markers)

# get celltypes
celltype_anno_cluster = AgrpPomcNTS_processed@meta.data %>% dplyr::group_by(preliminary_clusters,celltype_annotation) %>%
  dplyr::count()
top_2_markers = top_2_markers %>% dplyr::left_join(celltype_anno_cluster,by=c("cluster"="preliminary_clusters"))
top_2_markers$cluster_annotation[top_2_markers$celltype_annotation != "Neurons"] =  paste0(top_2_markers$cluster[top_2_markers$celltype_annotation != "Neurons"],": ",top_2_markers$celltype_annotation[top_2_markers$celltype_annotation != "Neurons"])
## add
temp_meta = dplyr::left_join(AgrpPomcNTS_processed@meta.data,top_2_markers[,c("cluster","cluster_annotation")],by=c("preliminary_clusters"="cluster"))
rownames(temp_meta) = temp_meta$Cell_ID

AgrpPomcNTS_processed@meta.data = temp_meta

##########
### DEGs - loaded from step 3
##########

DEG_ConditionAgrp = data.table::fread(paste0(project_path,"AgrpPomcNTS_DEG_Condition2_WT_CNO.txt"),data.table = F)
total_genes = length(unique(DEG_ConditionAgrp$gene))
DEG_ConditionAgrp$pval = DEG_ConditionAgrp[,grepl("p_Cond",colnames(DEG_ConditionAgrp))]
DEG_ConditionAgrp = DEG_ConditionAgrp %>% dplyr::group_by(cluster) %>% dplyr::mutate(pval_adj = p.adjust(pval,method="holm")) %>% dplyr::ungroup()
DEG_ConditionAgrp$cluster = as.character(DEG_ConditionAgrp$cluster)

DEG_ConditionAgrp_filt = DEG_ConditionAgrp[abs(DEG_ConditionAgrp$avg_log2FC) > 0.1 & DEG_ConditionAgrp$pval < 0.1, ]

data.table::fwrite(DEG_ConditionAgrp_filt,paste0(project_path,"deg_nebula_all.txt"),sep="\t")

##########
### Summarize results
##########

top_iegs = c("Fos","Fosl2","Homer1","Nr4a3","Nr4a1","Gem","Jun","Junb","Btg1","1700016P03Rik")
fishersMethod = function(x){pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)}

# The pvalues of 11 relevant IEGs were combined for each cluster, only including IEGs expressed in at least 5% of cells and significantly different (P<0.05) in the cluster), using Fisher's method (https://doi.org/10.2307%2F2681650).
# We then used the the negative logarithm of the chi-squared test statistic as a score to indicate neuronal activation. If the average fold changes of significant IEGs was negative we changed the sign of the score, thereby indicating potential inihibition.

# subset DEGs to IEGs
top_iegs_Activation = DEG_ConditionAgrp %>%
  dplyr::filter(gene %in% top_iegs & pval < 0.05) %>%
  dplyr::mutate(log10_pval = -log10(pval)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::mutate(fisher_pval = fishersMethod(pval),fisher_score =-log(fisher_pval)) %>%
  dplyr::mutate(mean_fc = mean(avg_log2FC)) %>%
  dplyr::arrange(desc(fisher_score))

## with correction:
# top_iegs_Activation = DEG_ConditionAgrp %>% 
#   dplyr::filter(gene %in% top_iegs ) %>%
#   dplyr::group_by(cluster) %>%
#   dplyr::add_count(name="n_tested_iegs" ) %>%
#   dplyr::mutate(pval_adj_ieg = p.adjust(pval,method = "holm",n = n_tested_iegs[1])) %>%
#   dplyr::mutate(log10_pval = -log10(pval)) %>% 
#   dplyr::filter(pval_adj_ieg < 1) %>%
#   dplyr::mutate(fisher_pval = fishersMethod(pval_adj_ieg),fisher_score =-log(fisher_pval)) %>%
#   dplyr::mutate(mean_fc = mean(avg_log2FC)) %>%
#   dplyr::arrange(desc(fisher_score))

top_iegs_Activation$cluster = as.character(top_iegs_Activation$cluster)

# get significant IEGs
top_iegs_Activation_sig = top_iegs_Activation %>%
  dplyr::filter(fisher_score > -log(0.05)) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(IEGs = paste0(gene,collapse = "/"))
top_iegs_Activation_sig$cluster = as.character(top_iegs_Activation_sig$cluster)

# get annotations
cluster_anno = AgrpPomcNTS_processed@meta.data %>%
  dplyr::group_by(preliminary_clusters,celltype_annotation,cluster_annotation) %>%
  dplyr::count(name="ncells")
cluster_anno$preliminary_clusters = as.character(cluster_anno$preliminary_clusters)

# aclaulte IEG sum to rank activated cell types
per_cluster_sum = top_iegs_Activation %>%
  dplyr::filter(fisher_score > -log(0.05)) %>%
  dplyr::mutate(fisher_score = case_when( mean_fc >=0 ~ fisher_score,
                                          mean_fc < 0 ~ fisher_score * -1)) %>%
  dplyr::filter(abs(mean_fc) >= 0.3 ) %>%
  dplyr::distinct(fisher_score,cluster) %>% dplyr::rename(ieg_score = fisher_score) %>%
  dplyr::mutate(ieg_label = cluster) %>%
  dplyr::mutate(activation_category =
                  case_when(
                    abs(ieg_score) >=0 & ieg_score < 3*-log(0.05) ~ "weak activation",
                    abs(ieg_score) >=3*-log(0.05) & ieg_score < 3*-log(0.001)  ~ "medium activation",
                    abs(ieg_score) >= 3*-log(0.001) ~ "strong activation"
                  ))
per_cluster_sum$cluster = as.character(per_cluster_sum$cluster)
per_cluster_sum$activation_category[per_cluster_sum$ieg_score < 0] = gsub("activation","inhibition",per_cluster_sum$activation_category[per_cluster_sum$ieg_score < 0])


## join to activation_result
activation_result = dplyr::left_join(per_cluster_sum,top_iegs_Activation_sig,by=c("cluster"="cluster")) %>%
  dplyr::left_join(cluster_anno,by=c("cluster"="preliminary_clusters"))

## clean
activation_result_short = activation_result %>%
  dplyr::select(cluster_annotation,cluster,ieg_score,ieg_category=activation_category,IEGs,ncells) %>%
  dplyr::arrange(desc(ieg_score))

data.table::fwrite(activation_result_short,file = paste0(project_path,"AgrpPomcNTS_activatedClusters.txt"),sep="\t")

##########
### add to metadata
##########

AgrpPomcNTS_withIEG = AgrpPomcNTS_processed

temp_meta = AgrpPomcNTS_withIEG@meta.data %>% dplyr::left_join(activation_result_short[,c("cluster","ieg_score","ieg_category")],by=c("preliminary_clusters"="cluster")) %>% as.data.frame()
rownames(temp_meta) = temp_meta$Cell_ID

AgrpPomcNTS_withIEG@meta.data = temp_meta

##########
### save object for use inf figure creation script
##########

saveRDS(AgrpPomcNTS_withIEG,paste0(project_path,"AgrpPomcNTS_withIEG.rds"))
