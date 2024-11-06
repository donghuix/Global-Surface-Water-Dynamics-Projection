#!/bin/csh

#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --job-name=extract        ## job_name
#SBATCH --account=m3780           ## project_name 
#SBATCH --time=10:00:00           ## time_limit
#SBATCH --nodes=1                 ## number_of_nodes
#SBATCH --tasks-per-node=64       ## number_of_cores                                                                                              
#SBATCH --output=mat.stdout1      ## job_output_filename
#SBATCH --error=mat.stderr1       ## job_errors_filename

ulimit -s unlimited
module load matlab

matlab  -nodisplay -nosplash <Scripts_08_Read_Projections_cal12.m> Scripts_08_Read_Projections_cal12.log
