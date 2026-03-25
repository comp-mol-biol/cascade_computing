import flow
import numpy as np
import signac
import papermill as pm
import os
import glob
import pandas as pd
import re
import shutil
import MDAnalysis as mda
from MDAnalysis import transformations
import pickle
import math
import subprocess

from ..utils import toolbox_plots as ptools
from ..utils import toolbox_pp as pptools
from ..utils import toolbox_viz as vztools
from ..utils import toolbox_sim as simtools

from ..config import PATH_GIT, EMAIL, PATH_FUNCTIONS, PATH_TOOLS, NUM_WORKERS


class MyProject(flow.FlowProject):
    pass

# operation is only executed if all preconditions are met, and at at least one postcondition is not
# fake pre & post conditions:

@MyProject.label
def breaker(job):
    """
    A persistent false label to prevent operations from being marked as finished.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle.

    Returns
    -------
    bool
        Always returns False.
    """
    return False

@MyProject.label
def encourager(job):
    """
    A persistent true label for flow control.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle.

    Returns
    -------
    bool
        Always returns True.
    """
    return True


@MyProject.post(breaker)
@MyProject.operation
def post_processing(job):
    """
    Sets up the directory structure and saves essential topology and trajectory 
    metadata for post-processing.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    # paths
    output_path = job.document['output_path']
    print(output_path)
    label_new = str(job.document['params']['run'])
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
    os.makedirs(bookkeeping_path, exist_ok=True)
    path_script_new = job.path
    
    # load topology 
    top = job.document['top']
    top_mix = job.document['top_mix'] 

    # load topology with orginal indexing
    top_org = job.document['top_org']

    # and short intial trajectory
    xtcs = job.document['xtc']
    exmpl_traj = xtcs[0]

    # load domains
    df_domains = job.data['df_domains']
    df_sys_domains = job.data['df_sys_domains']

    job.doc.domains = pptools.find_and_sort_min_max_columns(job.data['df_domains'])

    # save xtcs file list
    with open(bookkeeping_path + '/' + label_new + '_xtc_list.pkl', 'wb') as file:
        pickle.dump(xtcs, file)
    
    # save post_processing
    shutil.copy(top, bookkeeping_path + '/' + label_new + "_md_run.tpr")
    shutil.copy(top_mix, bookkeeping_path + '/' + label_new + "_topol_mix.top")

    # save molecular domains
    df_domains.to_parquet(bookkeeping_path + '/' + label_new + '_domains.parquet')
    df_sys_domains.to_parquet(bookkeeping_path + '/' + label_new + '_sys_domains.parquet')


@MyProject.post(breaker)
@MyProject.operation
def transform(job):
    """
    Prepares and submits a SLURM batch script to run the trans-pi-analysis 
    Python script.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    output_path = job.document['output_path']
    label_new = str(job.document['params']['run'])
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
    os.makedirs(bookkeeping_path, exist_ok=True)
    path_script_new = job.path

    # python script
    str_run_contacts = "python3 -u -m " + "src.compute.functions.traj_transform_pi --output " + output_path + " --label " + label_new + " --step " + str(1) #+ " --frames " + str(0)
    
    print(str_run_contacts)
    l_commands = ['export PYTHONPATH=$PYTHONPATH:'+PATH_GIT,'cd '+PATH_GIT,str_run_contacts]

    content = pptools.create_bash_file(l_commands, hours=1, nodes=1, name="trans_" + label_new)

    file_trans = path_script_new + '/' + label_new + "_run_trans.sh"
    log_trans = path_script_new + '/' + '%j_' + label_new + '_trans.out'
    with open(file_trans, 'w') as file:
            file.write(content)
    str_run = "! sbatch --mail-type=ALL --mail-user=" + EMAIL + " --output=" + log_trans + ' ' + file_trans
    os.system(str_run)


@MyProject.post(breaker)
@MyProject.operation
def contacts(job):
    """
    Defines contact calculation parameters, creates a job metadata parquet file,
    and submits a high-memory batch job for contact analysis.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    label_new = str(job.document['params']['run'])
    path_script_new = job.path
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
  
    job.document['cont_params'] = {
        "d_step": 1, 
        "delta_ha": 6, 
        "buffer": 0.4 * np.sqrt(3),
        "eps": 2,
        "bool_breaks": True,
        "breaks_tol": 1,
        "n_workers": 8, 
        "split_parts": 10,
        "traj_min": 0,
        "traj_max": 0,
        "bool_sc": True,
        "bool_bb": False,
        "pwi": False,
        "bool_q": False,
        "bool_bonds": False,
        "bool_debug": False
    }

    #c0=10
    #ck=20
    
    cont_params = job.document['cont_params']

    d_step      = cont_params["d_step"]
    delta_ha    = cont_params["delta_ha"]
    buffer      = cont_params["buffer"]
    eps         = cont_params["eps"]
    bool_breaks = cont_params["bool_breaks"]
    breaks_tol  = cont_params["breaks_tol"]
    n_workers   = cont_params["n_workers"]
    split_parts = cont_params["split_parts"]
    traj_min    = cont_params["traj_min"]
    traj_max    = cont_params["traj_max"]
    bool_sc     = cont_params["bool_sc"]
    bool_bb     = cont_params["bool_bb"]
    bool_pwi    = cont_params["pwi"]
    bool_q      = cont_params["bool_q"]
    bool_bonds  = cont_params["bool_bonds"]
    bool_debug  = cont_params["bool_debug"]

    # update metadata
    name = "minimal_bondss_ts_" + str(traj_min) + "_" + str(traj_max)
    job.document['contact_meta'] = name + '_accont_s' + str(d_step) + "_dha_" + str(delta_ha).replace('.', '_') + "_eps_" + str(eps) + "_brks_" + str(bool_breaks) + "_tol_" + str(breaks_tol) + "_buf" + str(buffer) + "_comp_" + str(n_workers) + "_" + str(split_parts) + "_sc_" + str(bool_sc) + "_bb_" + str(bool_bb) + "_pwi_" + str(bool_pwi) + "_bonds_" + str(bool_bonds) + "_debug_" + str(bool_debug).replace('.', '_')
    contact_path = bookkeeping_path + "/" + job.document['contact_meta']

    print("contact_path", contact_path)
    job.document['contact_path'] = contact_path
    os.makedirs(contact_path, exist_ok=True)

    job_dict = {'sgnc_path': job.doc.params.sgnc_path, 'label': job.doc.params.run, 'contact_path': job.doc.contact_path, 'output_path': job.doc.output_path, "id": job.id}
    job_df = pd.DataFrame([job_dict])
    job_df.to_parquet(bookkeeping_path + "/" + label_new + "_job_dict.parquet", index=False)
    
    # create bash file
    log_cont = contact_path + '/' + '%j_con_' + label_new + '.out'
    str_run_contacts = "python3 -u -m " + "src.compute.functions.calc_contacts_opt --path "+ bookkeeping_path+" --label "+label_new+" --output "+contact_path+" --cutoff_ha "+str(delta_ha)+" --step "+str(d_step)+" --cutoff_mol "+str(buffer)+" --cutoff_eps "+ str(eps)+" --breaks "+str(bool_breaks)+" --breaks_tol "+str(breaks_tol)+" --n_workers "+str(n_workers)+" --split_parts "+str(split_parts)+" --traj_min "+str(traj_min)+" --traj_max "+ str(traj_max)+" --sidechains "+ str(bool_sc)+" --backbone "+ str(bool_bb)+" --pwi "+ str(bool_pwi)+" --q "+ str(bool_q)+" --bonds "+ str(bool_bonds)+" --debug "+str(bool_debug)#+" --c0 "+ str(c0)+" --ck "+ str(ck)
    print(str_run_contacts)
    l_commands = [
        'export PYTHONPATH=$PYTHONPATH:'+PATH_GIT,
        'cd '+PATH_GIT,
        'python -c "import bokeh; print(\'Bokeh version:\', bokeh.__version__)"',
        'python -c "import pandas; print(\'Pandas version:\', pandas.__version__)"',
        'python -c "import dask, distributed; print(\'Dask:\', dask.__version__, \'| Distributed:\', distributed.__version__)"',
        str_run_contacts
    ]

    run_cont_file = contact_path + '/' + label_new + "_calc_contacts.sh"
    content = pptools.create_bash_file_memory(l_commands, hours=15, nodes=1, name="cc_contacts_" + label_new)
    #content = pptools.create_bash_file_high_memory(l_commands, hours=1, nodes=1, name="cc_contacts_" + label_new)
    with open(run_cont_file, 'w') as file:
        file.write(content)

    # run
    str_con = '! sbatch --mail-type=ALL --mail-user=' + EMAIL + ' --output=' + log_cont + " " + run_cont_file
    os.system(str_con)



@MyProject.post(breaker)
@MyProject.operation
def contacts_splid(job):
    """
    Defines contact calculation parameters, creates a job metadata parquet file,
    and submits a high-memory batch job for contact analysis.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    label_new = str(job.document['params']['run'])
    path_script_new = job.path
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
  
    job.document['cont_params'] = {
        "d_step": 1, #10 --> 1ns in MRT
        "delta_ha": 6, 
        "buffer": 0.4 * np.sqrt(3),
        "eps": 2,
        "bool_breaks": True,
        "breaks_tol": 1,
        "n_workers": 20, #40
        "split_parts": 40,
        "traj_min": 0,
        "traj_max": 0,
        "bool_sc": True,
        "bool_bb": False,
        "pwi":False,
        "bool_q":False,
        "bool_bonds": True,
        "bool_debug":True
    }

    cont_params = job.document['cont_params']

    d_step      = cont_params["d_step"]
    delta_ha    = cont_params["delta_ha"]
    buffer      = cont_params["buffer"]
    eps         = cont_params["eps"]
    bool_breaks = cont_params["bool_breaks"]
    breaks_tol  = cont_params["breaks_tol"]
    n_workers   = cont_params["n_workers"]
    split_parts = cont_params["split_parts"]
    traj_min    = cont_params["traj_min"]
    traj_max    = cont_params["traj_max"]
    bool_sc     = cont_params["bool_sc"]
    bool_bb     = cont_params["bool_bb"]
    bool_pwi    = cont_params["pwi"]
    bool_q      = cont_params["bool_q"]
    bool_bonds  = cont_params["bool_bonds"]
    bool_debug  = cont_params["bool_debug"]

    # update metadata
    name = "contacts_type_ts_" + str(traj_min) + "_" + str(traj_max)
    job.document['contact_meta'] = name + '_accont_s' + str(d_step) + "_dha_" + str(delta_ha).replace('.', '_') + "_eps_" + str(eps) + "_brks_" + str(bool_breaks) + "_tol_" + str(breaks_tol) + "_buf" + str(buffer) + "_comp_" + str(n_workers) + "_" + str(split_parts) + "_sc_" + str(bool_sc) + "_bb_" + str(bool_bb) + "_pwi_" + str(bool_pwi) + "_bonds_" + str(bool_bonds) + "_debug_" + str(bool_debug).replace('.', '_')
    contact_path = bookkeeping_path + "/" + job.document['contact_meta']

    print("contact_path", contact_path)
    job.document['contact_path'] = contact_path
    os.makedirs(contact_path, exist_ok=True)

    job_dict = {'sgnc_path': job.doc.params.sgnc_path, 'label': job.doc.params.run, 'contact_path': job.doc.contact_path, 'output_path': job.doc.output_path, "id": job.id}
    job_df = pd.DataFrame([job_dict])
    job_df.to_parquet(bookkeeping_path + "/" + label_new + "_job_dict.parquet", index=False)
    
    
    l_splid_parts=[[10,14], [14,17], [17,20]]

    for c0, ck in l_splid_parts:

        # run command
        str_run_contacts = "python3 -u " + PATH_FUNCTIONS + "/calc_contacts_opt.py --path "+ bookkeeping_path+" --label "+label_new+" --output "+contact_path+" --cutoff_ha "+str(delta_ha)+" --step "+str(d_step)+" --cutoff_mol "+str(buffer)+" --cutoff_eps "+ str(eps)+" --breaks "+str(bool_breaks)+" --breaks_tol "+str(breaks_tol)+" --n_workers "+str(n_workers)+" --split_parts "+str(split_parts)+" --traj_min "+str(traj_min)+" --traj_max "+ str(traj_max)+" --sidechains "+ str(bool_sc)+" --backbone "+ str(bool_bb)+" --pwi "+ str(bool_pwi)+" --q "+ str(bool_q)+" --bonds "+ str(bool_bonds)+" --debug "+str(bool_debug)+" --c0 "+ str(c0)+" --ck "+ str(ck)
        
        print(str_run_contacts)
        l_commands = [
            'python -c "import bokeh; print(\'Bokeh version:\', bokeh.__version__)"',
            'python -c "import pandas; print(\'Pandas version:\', pandas.__version__)"',
            'python -c "import dask, distributed; print(\'Dask:\', dask.__version__, \'| Distributed:\', distributed.__version__)"',
            str_run_contacts
        ]

        log_cont=contact_path +'/'+ '%j_con_'+ str(c0)+"_"+str(ck)+"_"+label_new+'.out'
        run_cont_file=contact_path+'/'+label_new+"_" + str(c0)+"_"+ str(ck) +"_calc_contacts.sh"
        nametag=str(c0)+"_"+ str(ck)
        content = pptools.create_bash_file_high_memory(l_commands, hours=5, nodes=1, name=nametag+"_cc_contacts_" + label_new)
        with open(run_cont_file, 'w') as file:
            file.write(content)

        # run
        str_con = '! sbatch --mail-type=ALL --mail-user=' + EMAIL + ' --output=' + log_cont + " " + run_cont_file
        os.system(str_con)

@MyProject.operation
def vizualization(job):
    """
    Generates visualization files as VMD TCL scripts for 
    system-wide and molecule-specific domain visualizations.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    label_new = str(job.document['params']['run'])
    output_path = str(job.document['output_path']) + '/postprocessing'
    
    # Clean up old artifacts
    for pattern in ['*.mp4', '*.gif']:
        files = glob.glob(os.path.join(output_path, pattern))
        for fi in files:
            try:
                os.remove(fi)
            except:
                pass

    try:
        domains_sys = job.data['df_sys_domains']
        
        content = vztools.create_vmd_file(domains_sys, './' + label_new + '_pi_clean.gro', 
                                          traj='./' + label_new + '_full_pi_ref.xtc', 
                                          lab=label_new + "_dom", type='VDW', 
                                          step=100, bool_dom=True, l_domains=job.doc.domains)

        print(content)
        with open(output_path + '/' + label_new + "_vmd_sys.tcl", 'w') as file:
            file.write(content)
        print(output_path + '/' + label_new + "_vmd_sys.tcl")
        content_t = vztools.create_vmd_file(domains_sys, './' + label_new + '_pi_clean.gro', 
                                            traj='./' + label_new + '_full_pi_ref.xtc', 
                                            lab=label_new + "_mol", type='VDW', 
                                            step=100, bool_dom=False, l_domains=job.doc.domains)

        print(content_t)

        
        with open(output_path + '/' + label_new + "_vmd_sys_t.tcl", 'w') as file:
            file.write(content_t) 

        print(output_path + '/' + label_new + "_vmd_sys_t.tcl")
    
        log_vcont = output_path + '/' + '%j_' + label_new + '_viz.out'
        str_run = "! sbatch --mail-type=ALL --mail-user=" + EMAIL + " --output=" + log_vcont + ' ' + PATH_FUNCTIONS + '/contact_viz.sh ' + output_path + " " + output_path + '/' + label_new + "_vmd_sys.tcl"
        
        log_vcont2 = output_path + '/' + '%j_' + label_new + '_viz2.out'
        str_run2 = "! sbatch --mail-type=ALL --mail-user=" + EMAIL + " --output=" + log_vcont2 + ' ' + PATH_FUNCTIONS + '/contact_viz.sh ' + output_path + " " + output_path + '/' + label_new + "_vmd_sys_t.tcl"
    
        os.system(str_run)
        os.system(str_run2)
    except:
        pass


@MyProject.post(breaker)
@MyProject.operation
def eval_contacts(job):
    """
    Evaluates computed contacts by running a contact evaluation script and 
    storing the resulting metadata in a job dictionary parquet file.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """

    bool_bonds=False
    n_workers=NUM_WORKERS-2

    label_new = str(job.document['params']['run'])
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
    contact_path = job.document['contact_path']
    result_folder = 'result_files'
    job.document['results_path'] = job.document['contact_path'] + '/' + result_folder

    job_dict = {'sgnc_path': job.doc.params.sgnc_path, 'label': job.doc.params.run, 'contact_path': job.doc.contact_path, 'results_path': job.doc.results_path, 'output_path': job.doc.output_path, "id": job.id, 'domains_ordered': ";".join(job.doc.domains)}
    job_df = pd.DataFrame([job_dict])
    job_df.to_parquet(bookkeeping_path + "/" + label_new + "_job_dict.parquet", index=False)
  
    log_cont = contact_path + '/' + '%j_con_eval_' + label_new + '.out'
    str_eval_contacts = "python3 -m " + "src.compute.functions.eval_contacts" + " --path " + bookkeeping_path + " --label " + label_new + " --output " + contact_path + " --result " + result_folder + " --pwi " + str(False) + " --path_tools " + PATH_TOOLS + " --bonds " + str(bool_bonds) + " --n_workers " + str(n_workers) 

    l_commands = [
        'export PYTHONPATH=$PYTHONPATH:'+PATH_GIT,
        'cd '+PATH_GIT,
        str_eval_contacts]
    run_cont_eval_file = contact_path + '/' + label_new + "_eval_contacts.sh"
    #content = pptools.create_bash_file_high_memory(l_commands, hours=1, nodes=1, name="eval_cc_" + label_new)
    content = pptools.create_bash_file_memory(l_commands, hours=1, nodes=1, name="eval_cc_" + label_new)
    with open(run_cont_eval_file, 'w') as file:
        file.write(content)

    str_con = '! sbatch --mail-type=ALL --mail-user=' + EMAIL + ' --output=' + log_cont + " " + run_cont_eval_file
    os.system(str_con)


@MyProject.post(breaker)
@MyProject.operation
def analysis(job):
    """
    Triggers the contact analysis script to process results generated in 
    the evaluation step.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    
    bool_bonds=False
    bool_pwi=False

    label_new = str(job.document['params']['run'])
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
    contact_path = job.document['contact_path']
    result_folder = 'result_files'

    str_ana_contacts = "python3 -m " + "src.compute.functions.pivot_analysis" + " --path_tools " + PATH_TOOLS + " --path " + bookkeeping_path + " --label " + label_new + " --output " + contact_path + " --result " + result_folder + " --pwi " + str(bool_pwi)  + " --bonds " + str(bool_bonds) 
    l_commands = [ 'export PYTHONPATH=$PYTHONPATH:'+PATH_GIT,
        'cd '+PATH_GIT,str_ana_contacts]

    run_cont_ana_file = contact_path + '/' + label_new + "_ana_contacts.sh"
    content = pptools.create_bash_file(l_commands, hours=1, nodes=1, name="ana_cc_" + label_new)
    
    with open(run_cont_ana_file, 'w') as file:
        file.write(content)

    log_cont = contact_path + '/' + '%j_con_ana_' + label_new + '.out'
    str_con = '! sbatch --mail-type=ALL --mail-user=' + EMAIL + ' --output=' + log_cont + " " + run_cont_ana_file
    os.system(str_con)


@MyProject.post(breaker)
@MyProject.operation
def derive_pivot(job):
    """
    Executes a pivot-specific evaluation script for the contact data.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    result_folder = 'result_files'
    label_new = str(job.document['params']['run'])
    bookkeeping_path = job.document['output_path'] + '/postprocessing'
    contact_path = job.document['contact_path']
    job.document['results_path'] = job.document['contact_path'] + '/' + result_folder

    job_dict = {'sgnc_path': job.doc.params.sgnc_path, 'label': job.doc.params.run, 'contact_path': job.doc.contact_path, 'results_path': job.doc.results_path, 'output_path': job.doc.output_path, "id": job.id, 'domains_ordered': ";".join(job.doc.domains)}
    job_df = pd.DataFrame([job_dict])
    job_df.to_parquet(bookkeeping_path + "/" + label_new + "_job_dict.parquet", index=False)
  
    log_cont = contact_path + '/' + '%j_con_eval_' + label_new + '.out'
    str_eval_contacts = "python3 " + PATH_FUNCTIONS + "/eval_pivots_only.py" + " --path " + bookkeeping_path + " --label " + label_new + " --output " + contact_path + " --result " + result_folder + " --pwi " + str(False)

    l_commands = [str_eval_contacts]
    run_cont_eval_file = contact_path + '/' + label_new + "_eval_contacts.sh"
    content = pptools.create_bash_file(l_commands, hours=1, nodes=1, name="tbc_eval_cc_" + label_new)
    with open(run_cont_eval_file, 'w') as file:
        file.write(content)

    str_con = '! sbatch --mail-type=ALL --mail-user=' + EMAIL + ' --output=' + log_cont + " " + run_cont_eval_file
    os.system(str_con)


@MyProject.post(breaker)
@MyProject.operation
def archiving(job):
    """
    Compresses the output folder into a .tar.gz archive using pigz.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle 
    """
    try:
        output_path = str(job.document['output_path'])
        output_folder = os.path.basename(output_path)
        archiv_folder = os.path.dirname(os.path.dirname(output_path))
        
        file_arch = os.path.join(archiv_folder, f"{job.document['params']['run']}.tar.gz")
        com = f"tar cf - -C {os.path.dirname(output_path)} {output_folder} | pigz -p 8 > {file_arch}"
        
        content = pptools.create_bash_file(com, hours=10, nodes=1, name='a_' + job.id)
        file_vcont = job.document['output_path'] + '/postprocessing/' + str(job.document['params']['run']) + "_run_acv.sh"
        log_vcont = job.document['output_path'] + '/postprocessing/' + '%j_' + str(job.document['params']['run']) + '_run_acv.out'
        with open(file_vcont, 'w') as file:
                file.write(content)
        str_run = "! sbatch --mail-type=ALL --mail-user=" + EMAIL + " --output=" + log_vcont + ' ' + file_vcont
        os.system(str_run)
    except:
        pass


@MyProject.post(breaker)
@MyProject.operation
def collecting_data(job):
    """
    Updates the job document with definitive paths to output files and contact results.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle 
    """
    prefix = job.document['output_path'] + '/postprocessing/'+ str(job.document['params']['run'])
    res_prefix = job.document['results_path'] + '/' + str(job.document['params']['run'])

    job.document['gro'] = prefix + '_pi_clean.gro'
    job.document['pdb'] = prefix + '_pi_clean.pdb'
    job.document['xtc'] = prefix + '_full_pi_ref.xtc'
    job.document['tcl'] = prefix + '_vmd_sys_t.tcl'

    job.document['contact_record'] = res_prefix + '_contacts_only_t.parquet'
    job.document['contact_meta'] = os.path.dirname(os.path.dirname(res_prefix))+  "/"+str(job.document['params']['run']) +'_contacts_res_meta.parquet'
    job.document['trajectory_meta'] =os.path.dirname( os.path.dirname(res_prefix))+"/"+ str(job.document['params']['run']) + '_time_status.parquet' 
    # Pivot files
    for res in ['prot', 'res_type', 'res_org', 'struc_id', 'res_dom']:
        job.document[f'{res}_prob'] = f'{res_prefix}_{res}_prob_pivot.parquet'
        job.document[f'{res}_freq'] = f'{res_prefix}_{res}_2time_count_pivot.parquet'
        job.document[f'{res}_perc'] = f'{res_prefix}_{res}_list_percistance2_distribution_pivot.parquet'

    exists1 = os.path.exists(job.document['gro'])
    exists2 = os.path.exists(job.document['pdb'])
    exists3 = os.path.exists(job.document['xtc'])
    exists4 = os.path.exists(job.document['tcl'])
    exists5 = os.path.exists(job.document['trajectory_meta'])
    
    print('gro',exists1,'pdb',exists2,'xtc', exists3, 'tcl',exists4, 'traj_meta',exists5 )


@MyProject.post(breaker)
@MyProject.operation
def connecting_data(job):
    """
    Updates the job document with definitive paths to output files and contact results.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle 
    """
    prefix = job.document['data_path'] + '/'+str(job.document['params']['run'])
    #res_prefix = job.document['results_path'] + '/' + str(job.document['params']['run'])
    print(prefix)

    job.document['gro'] = prefix + '_pi_clean.gro'
    job.document['pdb'] = prefix + '_pi_clean.pdb'
    job.document['xtc'] = prefix + '_full_pi_ref.xtc'
    job.document['tcl'] = prefix + '_vmd_sys_t.tcl'

    job.document['contact_record'] = prefix + '_contacts_only_t.parquet'
    job.document['contact_meta'] = prefix+'_contacts_res_meta.parquet'
    job.document['trajectory_meta'] =prefix+ '_time_status.parquet' 
    # Pivot files
    for res in ['prot', 'res_type', 'res_org', 'struc_id', 'res_dom']:
        job.document[f'{res}_prob'] = f'{prefix}_{res}_prob_pivot.parquet'
        job.document[f'{res}_freq'] = f'{prefix}_{res}_time_count_pivot.parquet'
        job.document[f'{res}_perc'] = f'{prefix}_{res}_list_percistance_distribution_pivot.parquet'

    exists1 = os.path.exists(job.document['gro'])
    exists2 = os.path.exists(job.document['pdb'])
    exists3 = os.path.exists(job.document['xtc'])
    exists4 = os.path.exists(job.document['tcl'])
    exists5 = os.path.exists(job.document['trajectory_meta'])
    
    print('gro',exists1,'pdb',exists2,'xtc', exists3, 'tcl',exists4, 'traj_meta',exists5 )



@MyProject.post(breaker)
@MyProject.operation
def document_files(job):
    
    target_folder =job.path+'/results_files_lifetimes'
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)
        print(f"Created directory: {target_folder}")
    
    # We put the paths into a list and copy them one by one
    files_to_copy = [
        job.document['gro'],
        job.document['pdb'],
        #job.document['xtc'],
        job.document['tcl'],
        job.document['trajectory_meta'],
        job.document['contact_record'],
        job.document['contact_meta']
    ]
    
    for res in ['res_type']:#['prot', 'res_type', 'res_org', 'struc_id', 'res_dom']:
        files_to_copy.append(job.document[f'{res}_prob'] )
        files_to_copy.append(job.document[f'{res}_freq'] )
        files_to_copy.append( job.document[f'{res}_perc'] )
    
    for file_path in files_to_copy:
        if os.path.exists(file_path):
            shutil.copy(file_path, target_folder)
            print(f"Successfully copied: {os.path.basename(file_path)}")
        else:
            print(f"WARNING: File not found: {file_path}")


@MyProject.post(breaker)
@MyProject.operation
def archiving_NOMAD(job):
    """
    takes the data from collecting_data and prepares 
    NOMAD entry submission.

    Parameters
    ----------
    job : signac.contrib.job.Job
        The job handle
    """
    print(job.document['cont_params'])
    print(job.id)
    print(job.sp)
    print(job.data['df_domains'])
    print(job.data['df_sys_domains'])
