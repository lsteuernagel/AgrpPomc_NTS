
# start
message(Sys.time(),": Load parameters and packages ")
message(".libPaths(): ")
message(paste0(.libPaths(),collapse = " | "))
library(scUtils,lib.loc = "/beegfs/scratch/bruening_scratch/lsteuernagel/R/user_lib/x86_64-pc-linux-gnu-library/4.2")

options(future.globals.maxSize= 10000*1024^2)

# get params-filename from commandline
command_args<-commandArgs(TRUE)
param_file = command_args[1]
# read all parameters and filepaths
parameter_list = jsonlite::read_json(param_file)
# if some fields are lists --> unlist
parameter_list = lapply(parameter_list,function(x){if(is.list(x)){return(unlist(x))}else{return(x)}})

# get path to functions
function_path = command_args[2]

# run function
message(Sys.time(),": Run differential_gene_expression ")
source(paste0(function_path,"R/pipeline_functions.R"))
source(paste0(function_path,"R/process_functions.R"))
deg_result = differential_gene_expression(param_file, verbose =TRUE)

# print(deg_result)
# saveRDS(deg_result,paste0(parameter_list$project_path,parameter_list$project_name,"DEG_",parameter_list$primary_variable,"_",parameter_list$reference_level,"_",parameter_list$res_name_add,".rds"))
# print("------2")

# save result
if(parameter_list$res_name_add == ""){
  filename = paste0(parameter_list$project_path,parameter_list$project_name,"_DEG_",parameter_list$primary_variable,"_",parameter_list$reference_level,".txt")
}else{
  filename = paste0(parameter_list$project_path,parameter_list$project_name,"_DEG_",parameter_list$primary_variable,"_",parameter_list$reference_level,"_",parameter_list$res_name_add,".txt")
}
data.table::fwrite(deg_result,filename,sep="\t")

# end
message(Sys.time(),": Complete ")

