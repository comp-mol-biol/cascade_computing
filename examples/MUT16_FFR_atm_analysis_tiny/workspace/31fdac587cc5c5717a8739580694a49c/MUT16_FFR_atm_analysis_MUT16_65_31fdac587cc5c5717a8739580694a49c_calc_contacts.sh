#!/bin/bash
#-----------------------------------------------------------------
# cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c
# 
#-----------------------------------------------------------------

#SBATCH -J "cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c"         	#todo Job name #todo
#SBATCH -o "cc_MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c".%j.out  	#todo Specify stdout output file (%j expands to jobId)
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
python3 -u /home/lubaltz/code/SFB1551/r08/GMX_models/py_tools/calc_contacts_2_11_tested_MG_run.py --path /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_10/dynamics/postprocessing --label MUT16_FFR_atm_analysis_MUT16_65_31fdac587cc5c5717a8739580694a49c --output /lustre/miifs01/project/m2_trr146/kugaurav/MUT16_FFR_atomistic/replica_10/dynamics/postprocessing/cf_pt_sc_ts_0_0_accont_s1_dha_3_8_eps_4.2_brks_True_tol_1_buf0.6928203230275509_comp_40_80_sc_True_bb_False_pwi_False_bonds_False --cutoff_ha 3.8 --step 1 --cutoff_mol 0.6928203230275509 --cutoff_eps 4.2 --breaks True --breaks_tol 1 --n_workers 40 --split_parts 80 --traj_min 0 --traj_max 0 --sidechains True --backbone False --pwi False --q False --bonds False


echo byebye bash
