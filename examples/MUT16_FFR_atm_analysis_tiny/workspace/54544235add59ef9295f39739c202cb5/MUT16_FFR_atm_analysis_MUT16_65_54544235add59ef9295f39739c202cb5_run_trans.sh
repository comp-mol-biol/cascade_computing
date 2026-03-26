#!/bin/bash
#-----------------------------------------------------------------
# ppMUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5
# 
#-----------------------------------------------------------------

#SBATCH -J "ppMUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5"         	#todo Job name #todo
#SBATCH -o "ppMUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel              	        # Partition/Queue name
#SBATCH -C skylake #broadwell #skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N 1                     	        #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t 20:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc             	   # Specify allocation to charge against
##SBATCH --mem=0
##SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
papermill -p output_path /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_4/dynamics /home/lubaltz/code/SFB1551/r08/GMX_models/analysis/t_trans-pi-analysis.ipynb /home/lubaltz/code/SFB1551/r08/projects/MUT16_FFR_atm_analysis/workspace/54544235add59ef9295f39739c202cb5/MUT16_FFR_atm_analysis_MUT16_65_54544235add59ef9295f39739c202cb5_trans-prot.ipynb


echo byebye bash
