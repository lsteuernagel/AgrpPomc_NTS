#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=240:00:00
#SBATCH --nice

# need to get 3 variables from call
# singularity image
singularity_image=$1
# $method_script : script that should be executed (depend on method)
method_script=$2
# $param_file : specifies which object, assay etc are used and which parameters of the method
param_file=$3
#function_path
function_path=$4

# Run script
echo "singularity exec "$singularity_image" Rscript "$method_script" "$param_file" "$function_path
srun singularity exec $singularity_image Rscript $method_script $param_file $function_path
