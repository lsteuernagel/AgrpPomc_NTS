
## 03: Run differential gene expression between conditions.

# Info: This step uses the slurm job version of the R scripts for semi automatic submission to the slurm queueing system
# All steps can all so be run directly in an active r session by using the functions in pipeline_functions.R ,e.g., scseq_pipeline and providing the parameters as a list. 

# need top level of pipeline
path_to_pipeline = "processing_pipeline/"

source(paste0(path_to_pipeline,"R/pipeline_functions.R"))
source(paste0(path_to_pipeline,"R/slurm_functions.R"))

## project settings
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"

# list scripts
list_slurm_scripts(slurm_dir = paste0(path_to_pipeline,"slurm/"))

condition_var = "Condition2"
param_file_name = paste0("agrpPomcNTS_deg_params_",condition_var,".json")
params_deg_Agrp = deg_parameters(param_file = param_file_name,
                                 project_path = project_path,
                                 input_seurat = paste0(project_path,"AgrpPomcNTS_processed.rds"),
                                 project_name = "AgrpPomcNTS",
                                 cluster_variable="preliminary_clusters",
                                 sample_variable="Sample_ID",
                                 primary_variable=condition_var,
                                 reference_level = "WT_CNO",
                                 other_variables=character(0),
                                 nCounts = "nCount_RNA",
                                 min_cells=10,
                                 min_pct_valid_samples = 0.5,
                                 min_pct = 0.05,
                                 assay="RNA",
                                 numCores = 56,
                                 padjust_method="bonferroni"
)
# submit to slurm
pipeline_id = submit_slurm(param_file_path = paste0(project_path,param_file_name),
                           script_path=paste0(path_to_pipeline,"slurm/slurm_differential_gene_expression.R"),
                           path_to_pipeline =  path_to_pipeline,
                           singularity_path = "~/Documents/r_scvi_v3_42.simg",
                           dependency_ids=NULL)