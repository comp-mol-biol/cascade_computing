#!/bin/bash
#-----------------------------------------------------------------
# eval_cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c
# 
#-----------------------------------------------------------------

#SBATCH -J "eval_cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c"         	#todo Job name #todo
#SBATCH -o "eval_cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel              	        # Partition/Queue name
#SBATCH -C broadwell #broadwell #skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N 1                     	        #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t 03:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc             	   # Specify allocation to charge against
#SBATCH --mem=100000
##SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
python3 /home/lubaltz/code/SFB1551/r08/GMX_models/py_tools/contact_eval_2_0_slim_multi_fixed_tested2_count_MG.py --path /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_10/dynamics/postprocessing --label MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c --output /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_10/dynamics/postprocessing/cf_ions_sc_ts_0_0_accont_s1_dha_5_eps_3_brks_True_tol_3_buf0.6928203230275509_comp_40_80_sc_True_bb_False_pwi_True_bonds_False --result result_files --pwi False


echo byebye bash
