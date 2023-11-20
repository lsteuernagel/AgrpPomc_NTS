
# start
message(Sys.time(),": Load parameters and packages ")

options(future.globals.maxSize= 10000*1024^2)
message(".libPaths(): ")
message(paste0(.libPaths(),collapse = " | "))
library(scUtils,lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.2")

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

if(is.null(parameter_list$res_name_add)){
  parameter_list$res_name_add = ""
}

# get path to functions
function_path = command_args[2]

# run function
message(Sys.time(),": Run scseq_pipeline ")
source(paste0(function_path,"R/pipeline_functions.R"))
source(paste0(function_path,"R/process_functions.R"))
processed_object = scseq_pipeline(param_file,verbose = TRUE,exec_singularity = FALSE )

# save result
saveRDS(processed_object,paste0(parameter_list$project_path,parameter_list$project_name,"_processed.rds"))

# end
message(Sys.time(),": Complete ")

