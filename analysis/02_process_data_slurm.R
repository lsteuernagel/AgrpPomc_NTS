
## 02: Run standard single cell QC and processing pipeline with scvi integration

# Info: This step uses the slurm job version of the R scripts for semi automatic submission to the slurm queueing system
# All steps can all so be run directly in an active r session by using the functions in pipeline_functions.R ,e.g., scseq_pipeline and providing the parameters as a list. 

# need top level of pipeline
path_to_pipeline = "processing_pipeline/" 

# source main functions
source(paste0(path_to_pipeline,"R/pipeline_functions.R"))
source(paste0(path_to_pipeline,"R/slurm_functions.R"))

## project settings -- adjust to local path !
raw_file_path = paste0("/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/","almu_2023_04_AgrpPomc_raw.rds")
param_file_name = "agrpPomcNTS_params.json"
project_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/data/2023-03-almu-snseq-analysis/AgrpPomc/processed/"

# list scripts
list_slurm_scripts(slurm_dir = paste0(path_to_pipeline,"slurm/"))

## pipeline parameters
params_scseq = scseq_parameters(param_file = param_file_name,
                                project_path = project_path,
                                input_seurat = raw_file_path,
                                project_name = "AgrpPomcNTS",
                                integration_name = "scvi_native",
                                batch_var = "Sample_ID",
                                minUMI = 800,
                                feature_set_size = 2000,
                                target_cluster_number = 250,
                                resolutions_to_try = c(0.5, 0.75, 1, 1.5, 2:20),
                                max_epochs = 80,
                                n_latent = 50,
                                doublet_rate = 0.02,
                                run_markers = TRUE)

# submit to slurm
pipeline_id = submit_slurm(param_file_path = paste0(project_path,param_file_name),
                           script_path=paste0(path_to_pipeline,"slurm/slurm_scseq_pipeline.R"),
                           path_to_pipeline = path_to_pipeline,
                           singularity_path = "~/Documents/r_scvi_v3_42.simg",
                           dependency_ids=NULL)
