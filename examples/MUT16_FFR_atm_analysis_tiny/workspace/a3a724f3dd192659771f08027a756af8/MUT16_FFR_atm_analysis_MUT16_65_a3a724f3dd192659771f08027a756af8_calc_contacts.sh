#!/bin/bash
#-----------------------------------------------------------------
# cc_MUT16_FFR_atm_analysis_MUT16_65_a3a724f3dd192659771f08027a756af8
# 
#-----------------------------------------------------------------

#SBATCH -J "cc_MUT16_FFR_atm_analysis_MUT16_65_a3a724f3dd192659771f08027a756af8"         	#todo Job name #todo
#SBATCH -o "cc_MUT16_FFR_atm_analysis_MUT16_65_a3a724f3dd192659771f08027a756af8".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p bigmem           # Partition/Queue name 'bigmem'
#SBATCH -C broadwell          # select either 'broadwell' or 'skylake'
#SBATCH -N 1           #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t 110:00:00     #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=0
#SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
python -c "import bokeh; print('Bokeh version:', bokeh.__version__)"
python -c "import pandas; print('Pandas version:', pandas.__version__)"
python -c "import dask, distributed; print('Dask:', dask.__version__, '| Distributed:', distributed.__version__)"
python3 -u /home/lubaltz/code/SFB1551/r08/GMX_models/py_tools/calc_contacts_2_11_tested_MG_run.py --path /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_8/dynamics/postprocessing --label MUT16_FFR_atm_analysis_MUT16_65_a3a724f3dd192659771f08027a756af8 --output /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_8/dynamics/postprocessing/cf_bb_ts_0_0_accont_s1_dha_6_eps_2_brks_True_tol_1_buf0.6928203230275509_comp_40_80_sc_False_bb_True_pwi_False_bonds_False --cutoff_ha 6 --step 1 --cutoff_mol 0.6928203230275509 --cutoff_eps 2 --breaks True --breaks_tol 1 --n_workers 40 --split_parts 80 --traj_min 0 --traj_max 0 --sidechains False --backbone True --pwi False --q False --bonds False


echo byebye bash
