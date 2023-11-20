
##########
### list_slurm_scripts
##########

#' List available slurm R scripts for execution
#'
#' All pipelien functions are wrapped in a slurm-ready Rscript that can be called from the sbatch script. This function lists all available scripts.
#' Names have to look like slurm_SCRIPTNAME.R
#'
#' @param slurm_dir where to look for scripts
#'
#' @return vector with script names
#'

list_slurm_scripts = function(slurm_dir = "slurm/"){
  all_files= list.files(path = slurm_dir,pattern = "slurm_.*\\.R",full.names = TRUE)
  return(gsub("//","/",all_files))
}


##########
### submit_slurm
##########

#' Submit Rscript job via sbatch to HPC
#'
#' This function submits a slurm job via system call
#' pipeline: Runs filtering, doublet detectiom, scvi, seurat, clustering, annotation and marker detection
#'
#' @param param_file_path full path to a json with all required parameters written with scseq_parameters
#' @param script_path path to slurm script
#' @param function_path where to find the pipeline
#' @param singularity_path singluarity image to use
#' @param jobname name for job. NULL with create automatically
#' @param dependency_ids slurm ids of jobs
#' @param outputfile name for log file. NULL with create automatically using log_path as dir
#' @param errorfile name for error log file. NULL with create automatically using log_path as dir
#' @param log_path if outputfile or errorfile are NULL, this path will be used for autmatic creation.
#' @param cpus_task how many cpus to allocate per node. NULL will default to mpi partition max
#' @param partition which mpi hpc partition to run on
#'
#' @return list of parameters for scseq_pipeline
#'

submit_slurm = function(
  param_file_path,
  script_path,
  path_to_pipeline = "",
  singularity_path = "~/Documents/r_scvi_v3_42.simg",
  jobname = NULL,
  dependency_ids = NULL,
  outputfile = NULL,
  errorfile = NULL,
  log_path = "/beegfs/scratch/bruening_scratch/lsteuernagel/slurm/other_logs/",
  cpus_task = NULL,
  partition = "blade-b"
) {

  #### Handle jobname
  if(is.null(jobname)){
    jobname = gsub("\\.R","",gsub("slurm_","",strsplit(x = script_path,split =  "/")[[1]][-1]))
  }

  #### Handle dependency_ids
  if(!is.null(dependency_ids)){
    dependencies = paste0(dependency_ids,collapse = ":")
  }else{
    dependencies = 0
  }
  #### Handle errorfile
  if(is.null(outputfile)){
    outputfile = paste0(log_path,"%j-script.out")
  }

  #### Handle errorfile
  if(is.null(errorfile)){
    errorfile = paste0(log_path,"%j-script.err")
  }

  #### Handle partition and cpus_task
  if(partition %in% c("blade-b","highmem")){
    if(is.null(cpus_task)){cpus_task = 56}else{cpus_task = min(cpus_task,56)}
  }else if(partition %in% c("blade-a")){
    if(is.null(cpus_task)){cpus_task = 40}else{cpus_task = min(cpus_task,40)}
  }else{
    stop("partition must be one of: blade-b, highmem, blade-a")
  }
  slurm_run_path = paste0(path_to_pipeline,"slurm/slurm_job.sh")
  function_path = path_to_pipeline
  #### Make sbatch call
  # output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile," --dependency=afterok:",dependencies,
  #                                " --kill-on-invalid-dep=yes ","--cpus-per-task ",cpus_task," --partition ",partition," slurm/slurm_job.sh ",
  #                                singularity_path," ",script_path," ",param_file_path),intern = TRUE)
  output_message = system(paste0("sbatch -J ",jobname," -o ",outputfile," -e ",errorfile,
                                 " --kill-on-invalid-dep=yes ","--cpus-per-task ",cpus_task," --partition ",partition," ",slurm_run_path," ",
                                 singularity_path," ",script_path," ",param_file_path," ",function_path),intern = TRUE)
  message("----- Slurm: ",output_message)
  slurm_id = gsub("Submitted batch job ","",output_message)

  #### return job id
  return(slurm_id)
}




