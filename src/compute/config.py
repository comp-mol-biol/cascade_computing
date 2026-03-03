#set PATHS
PATH_GIT ='/home/lubaltz/code/cascade_computing'
PATH_FUNCTIONS = PATH_GIT + '/src/compute/functions'
PATH_TOOLS = PATH_FUNCTIONS + '/src/compute/utils'

#set HPC variables
EMAIL = 'lubaltz@uni-mainz.de'
NUM_WORKERS = 64
OS_GB=50
LOCALSCRATCH="localscratch"

# parameters fraction of native contacts
BETA_CONST = 5  # 1/nm
LAMBDA_CONST = 1.2 # coarse grained 

#SLURM HEADER
slurm_header="""#!/bin/bash
#-----------------------------------------------------------------
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	
#SBATCH -o "{0}".%j.out  	
#SBATCH -p parallel              	        
#SBATCH -C broadwell     
#SBATCH -N {1:01}                     	       
#SBATCH -t {2:02}:00:00              	       
#SBATCH -A m2_komet331hpc            
#SBATCH --mem=246000
#SBATCH --exclusive

# Commands
"""

#ENVIROMENT
PATH_ENV='/home/lubaltz/p_3_10/bin/activate'