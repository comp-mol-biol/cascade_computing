import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.coordinates.XTC import XTCWriter
import pandas as pd
import os  
import glob
import numpy as np
import itertools
import argparse
from distutils.util import strtobool
import distributed
from distributed import Client, LocalCluster, get_worker
from dask import delayed
from dask.distributed import performance_report
import dask
import dask.dataframe as dd
import shutil
import psutil
import sys
import inspect
import re
import socket
import gc
import ctypes

# Path setup
path_tools='/home/lubaltz/code/cascade_computing' 
sys.path.append(path_tools + '/functions')
sys.path.append(path_tools + '/utils')
sys.path.append(path_tools + '/signac')
import toolbox_interactions as interactiontools


parser = argparse.ArgumentParser(description='Parallel Contact Map Calculation')
parser.add_argument('--path', action="store", dest='path', default='./')
parser.add_argument('--output', action="store", dest='output', default='./')
parser.add_argument('--label', action="store", dest='label', default='ohoh')
parser.add_argument('--cutoff_ha', action="store", dest='delta_ha', default=4.5)
parser.add_argument('--step', action="store", dest='d_step', default=1)
parser.add_argument('--cutoff_mol', action="store", dest='buffer', default=1)
parser.add_argument('--cutoff_eps', action="store", dest='delta_eps', default=1)
parser.add_argument('--breaks', action="store", dest='bool_breaks', default='True')
parser.add_argument('--breaks_tol', action="store", dest='eps_ts', default=1)
parser.add_argument('--n_workers', action="store", dest='n_workers', default=64)
parser.add_argument('--split_parts', action="store", dest='split_parts', default=10)
parser.add_argument('--traj_min', action="store", dest='traj_min', default=0)
parser.add_argument('--traj_max', action="store", dest='traj_max', default=0) 
parser.add_argument('--sidechains', action="store", dest='b_sidechains', default='False') 
parser.add_argument('--backbone', action="store", dest='b_bb', default='False') 
parser.add_argument('--pwi', action="store", dest='bool_pwi', default='False') 
parser.add_argument('--q', action="store", dest='bool_q', default='False') 
parser.add_argument('--bonds', action="store", dest='bool_bonds', default='False') 
parser.add_argument('--debug', action="store", dest='bool_debug', default='False')
parser.add_argument('--continue', action="store", dest='bool_filtering', default='False')
parser.add_argument('--c0', action="store", dest='c0', default=0)
parser.add_argument('--ck', action="store", dest='ck', default=0) 


# Parse Arguments
args = parser.parse_args()
path = args.path
output_path = args.output
label = args.label
delta_ha = float(args.delta_ha)  
delta_eps = float(args.delta_eps)   
d_step = int(args.d_step) 
buffer = float(args.buffer) 
bool_breaks = bool(strtobool(args.bool_breaks))
eps_ts = int(args.eps_ts)
n_workers = int(args.n_workers)
split_parts = int(args.split_parts)
traj_min = int(args.traj_min) 
traj_max = int(args.traj_max) 
b_sidechains = bool(strtobool(args.b_sidechains))
b_bb = bool(strtobool(args.b_bb))
bool_pwi = bool(strtobool(args.bool_pwi))
bool_q = bool(strtobool(args.bool_q))
bool_bonds = bool(strtobool(args.bool_bonds))
bool_debug = bool(strtobool(args.bool_debug))
bool_filtering = bool(strtobool(args.bool_filtering))
c0 = int(args.c0)
ck = int(args.ck)

if ck == 0 and c0 == 0:
    bool_splid = False
    print('full system - no splid')
else:
    bool_splid = True

# Global Parameters
BETA_CONST = 5  # 1/nm
LAMBDA_CONST = 1.2 # coarse grained 

# Set memory per worker
mem_str = str(psutil.virtual_memory().total / n_workers // 1024**3) + 'GB'
print('memory-per-worker', mem_str)


def cleanup_temp_folders(base_path, name='traj_parts', complete=False):
    """
    Remove temporary folders generated during the run.

    Parameters
    ----------
    base_path : str
        The root directory containing the temporary folders.
    name : str, optional
        Specific folder name to clean up. Default is 'traj_parts'.
    complete : bool, optional
        If True, removes all related temporary folders (tmp_parts, cont_map, etc.).
        If False, removes only the folder specified by `name`.
    """
    if complete:
        folder_patterns = [name, "tmp_parts", "cont_map_*", "lb_cont_map_*"]
    else:
        folder_patterns = [name]

    for pattern in folder_patterns:
        for folder in glob.glob(os.path.join(base_path, pattern)):
            if os.path.isdir(folder):
                try:
                    shutil.rmtree(folder, ignore_errors=True)
                    print(f"Deleted temporary folder: {folder}", flush=True)
                except Exception as e:
                    print(f"Could not delete {folder}: {e}", flush=True)


def split_trajectory_into_parts(top_file, traj_file, traj_min, traj_max, output_dir="split_trajectories", parts=10, d_step=1):
    """
    Split a large trajectory into smaller XTC chunks for parallel processing.

    Parameters
    ----------
    top_file : str
        Path to the topology file.
    traj_file : str
        Path to the trajectory file.
    traj_min : int
        Start frame index.
    traj_max : int
        End frame index.
    output_dir : str, optional
        Directory to save the chunks.
    parts : int, optional
        Number of chunks to create.
    d_step : int, optional
        Stride for skipping frames.

    Returns
    -------
    int
        The size (number of frames) of each chunk.
    """
    u = mda.Universe(top_file, traj_file, continuous=True, in_memory=False)
    print(len(u.trajectory), "length of input traj", flush=True)
    subset_trajectory = u.trajectory[traj_min:traj_max:d_step]
    n_frames = len(subset_trajectory)
    chunk_size = n_frames // parts

    os.makedirs(output_dir, exist_ok=True)
    print(f"Trajectory split into {parts} parts in '{output_dir} parts of total '{n_frames} frames with chunk size {chunk_size}'.")

    for i in range(parts):
        start = i * chunk_size
        end = (i + 1) * chunk_size if i < parts - 1 else n_frames
        print("from", start, end)
        output_file = os.path.join(output_dir, f"chunk_{i+1}.xtc")

        with XTCWriter(output_file, n_atoms=u.atoms.n_atoms, multiframe=True) as writer:
            for ts in subset_trajectory[start:end]: 
                writer.write(u.atoms)
        print(f"Written part {i+1}: frames {start} to {end - 1}")

    return chunk_size


def debug_partition_types(df, label):
    """
    Debug helper to inspect Dask partition data types and content.
    
    Parameters
    ----------
    df : pandas.DataFrame
    label : str
    
    Returns
    -------
    pandas.DataFrame
    """
    try:
        print(f"\n--- DEBUG [{label}] ---", flush=True)
        if len(df) > 0:
            sample_val = df['cont_a'].iloc[0]
            #print(f"Sample 'cont_a' type: {type(sample_val)}", flush=True)
            #print(f"Sample 'cont_a' content: {sample_val}", flush=True)

            sample_val2 = df['hbond'].iloc[0]
            #print(f"Sample 'hbond' type: {type(sample_val2)}", flush=True)
            #print(f"Sample 'hbond' content: {sample_val2}", flush=True)

            sample_val3 = df['salt_bridge'].iloc[0]
            #print(f"Sample 'salt_bridge' type: {type(sample_val3)}", flush=True)
            #print(f"Sample 'salt_bridge' content: {sample_val3}", flush=True)
            
            if isinstance(sample_val, str):
                print("!!! Data is a string here !!!", flush=True)
                if "..." in sample_val:
                    print("!!! ALERT: Data is truncated here !!!", flush=True)
            elif isinstance(sample_val, (np.ndarray, list)):
                print(">>> Data is a list/array here.", flush=True)
        else:
            print("Partition is empty.", flush=True)
    except Exception as e:
        print(f"Debug failed: {e}", flush=True)
    return df
        

def clean_and_parse(x):
    """
    Clean string representations of lists back into actual lists of integers.
    
    Parameters
    ----------
    x : str or list
    
    Returns
    -------
    list
    """
    if isinstance(x, list):
        print("yay list")
        return x
        
    if isinstance(x, str):
        if "..." in x:
            print(f"!!! TRUNCATION DETECTED !!!\nRaw Input: {x}", flush=True)
        
        tokens = x.strip("[]").replace(',', ' ').split()
        return [int(t) for t in tokens]
    return []


def clean_partition(df): 
    """
    Apply cleaning to contact columns in a dataframe partition.
    
    Parameters
    ----------
    df : pandas.DataFrame
    
    Returns
    -------
    pandas.DataFrame
    """
    df = df.copy()
    df['cont_a'] = df['cont_a'].apply(clean_and_parse)
    df['cont_b'] = df['cont_b'].apply(clean_and_parse)
    return df


def annotate_and_save_contacts(contacts_ddf, A_mol_sel, B_mol_sel, top, traj, a, b, path, name):
    """
    Annotate contact data with residue information and save to Parquet.
    
    Parameters
    ----------
    contacts_ddf : dask.dataframe.DataFrame
    A_mol_sel : str
    B_mol_sel : str
    top : str
    traj : str
    a : int
    b : int
    path : str
    name : str
    
    Returns
    -------
    tuple
        ([a, b], output_dir)
    """
    print('Im annoating here', flush=True)
    u = mda.Universe(top, traj, format='xtc', topology_format='pdb')
    A_mol = u.select_atoms(A_mol_sel)
    B_mol = u.select_atoms(B_mol_sel)

    # Build fast lookup dictionaries for residue ID and name
    id_to_res_a = {atom.id: res.resid for res in A_mol.residues for atom in res.atoms}
    id_to_type_a = {atom.id: res.resname for res in A_mol.residues for atom in res.atoms}
    
    id_to_res_b = {atom.id: res.resid for res in B_mol.residues for atom in res.atoms}
    id_to_type_b = {atom.id: res.resname for res in B_mol.residues for atom in res.atoms}

    # read contact-data
    meta_dict = {
        'frame': pd.Series([], dtype='int64'),
        'n_cont': pd.Series([], dtype='int64'),
        'cont_a': pd.Series([], dtype='object'),
        'cont_b': pd.Series([], dtype='object'),
        'q_values': pd.Series([], dtype='object'),
    }

    if bool_bonds:
        meta_dict.update({
            'cation_pi': pd.Series([], dtype='object'),
            'pi_stacking': pd.Series([], dtype='object'),
            'hbond': pd.Series([], dtype='object'),
            'salt_bridge': pd.Series([], dtype='object')
        })

    meta_contacts = pd.DataFrame(meta_dict)

    if bool_debug:
        contacts_ddf = contacts_ddf.map_partitions(debug_partition_types, label="Initial Load", meta=meta_contacts)
    else:
        contacts_ddf = contacts_ddf.map_partitions(lambda df: df, meta=meta_contacts)

    df = contacts_ddf[contacts_ddf['n_cont'] > 0].copy()

    meta_clean = df._meta.copy()
    meta_clean['cont_a'] = pd.Series(dtype=object)
    meta_clean['cont_b'] = pd.Series(dtype=object)
    
    df = df.map_partitions(clean_partition, meta=meta_clean)

    df = df.explode(['cont_a', 'cont_b'])
    df = df.repartition(npartitions=n_workers * 2)
    df = df.dropna(subset=['cont_a', 'cont_b'])

    df['cont_a'] = df['cont_a'].astype('int32')
    df['cont_b'] = df['cont_b'].astype('int32')

    df['res_a'] = np.nan
    df['type_a'] = np.nan
    df['frac_a'] = np.nan
    df['res_b'] = np.nan
    df['type_b'] = np.nan
    df['frac_b'] = np.nan
    
    # Vectorized mapping
    df['res_a'] = df['cont_a'].map(id_to_res_a, meta=('res_a', 'Int32'))
    df['type_a'] = df['cont_a'].map(id_to_type_a, meta=('type_a', 'object'))
    
    df['res_b'] = df['cont_b'].map(id_to_res_b, meta=('res_b', 'Int32'))
    df['type_b'] = df['cont_b'].map(id_to_type_b, meta=('type_b', 'object'))

    expected_cols = {'cont_a': 'int32', 'cont_b': 'int32'}
    for col, dtype in expected_cols.items():
        if col not in df.columns:
            df[col] = np.nan
            df[col] = df[col].astype(dtype)

    if len(df) == 0:
        print(f"Skipping file write for empty contact set: pair {a}-{b}")
        return ([a, b], None) 

    # Save
    out_dir = os.path.join(path, f"{name}_{a}_{b}")
    os.makedirs(out_dir, exist_ok=True)
    df.to_parquet(out_dir, engine="pyarrow", write_index=False)
  
    return ([a, b], out_dir)


def empty_contact_df():
    """Create an empty DataFrame with the correct schema for contacts."""
    meta_dict = {
        'frame': pd.Series([], dtype='int64'),
        'n_cont': pd.Series([], dtype='int64'),
        'cont_a': pd.Series([], dtype='object'),
        'cont_b': pd.Series([], dtype='object'),
        'q_values': pd.Series([], dtype='object'),
    }

    if bool_bonds:
        meta_dict.update({
            'cation_pi': pd.Series([], dtype='object'),
            'pi_stacking': pd.Series([], dtype='object'),
            'hbond': pd.Series([], dtype='object'),
            'salt_bridge': pd.Series([], dtype='object')
        })

    return pd.DataFrame(meta_dict)


def release_memory():
    """Manually trigger garbage collection and memory trim."""
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except Exception:
        pass


@delayed
def process_frame_chunk_for_split(top, traj, group_a_sel, group_b_sel, d_max, radius, k, chunk_size):
    """
    Process a chunk of frames to calculate contacts and specific interactions.
    
    Parameters
    ----------
    top, traj : str
    group_a_sel, group_b_sel : str
    d_max, radius : float
    k, chunk_size : int
    
    Returns
    -------
    pandas.DataFrame
    """
    if bool_debug:
        #check for correct file
        match = re.search(r'chunk_(\d+)', traj)
        if match:
            # Convert 1-based file index to 0-based calculation index
            k = int(match.group(1)) - 1
        else:
            print('check traj part', traj, flush=True)

        #check workload on worker
        worker_id = get_worker()
        print('part index',k, traj, flush=True)
        print('WORKER_ID', worker_id.name, flush=True)
        print(f"{worker_id.memory_manager.memory_limit / 1024**3:.2f} GB", flush=True)
        print(f"[{worker_id}] memory before processing {traj}: {psutil.Process().memory_info().rss / 1024**2:.2f} MB", flush=True)

    u = mda.Universe(top, traj, format='xtc', topology_format='pdb')

    str_a = group_a_sel + " and not name H*"
    u_group_a_sel = u.select_atoms(group_a_sel) 
    group_a = u_group_a_sel.select_atoms("not name H*")

    str_b = group_b_sel + " and not name H*"
    group_b = u.select_atoms(str_b) 
    
    if bool_bonds:
        group_a_sc_H = """
            ((resname SER and (name HG HG1)) or
            (resname THR and name HG1) or
            (resname TYR and name HH) or
            (resname ASN and (name HD21 HD22)) or
            (resname GLN and (name HE21 HE22)) or
            (resname HIS HSD HSE HSP and (name HD1 HE2)) or
            (resname TRP and name HE1) or
            (resname LYS and (name HZ1 HZ2 HZ3)) or
            (resname ARG and (name HE HH11 HH12 HH21 HH22)))
        """

        #for h-bonds - expand selection to hydrogens 
        #str_c=group_a_sel+ " and " +group_a_sc_H
        #group_c = u.select_atoms(str_c)
        #group_c=u_group_a_sel.select_atoms(group_a_sc_H)
        #print("c_selection", set(group_c.atoms.names), len(group_c.atoms.ids), flush=True)

        #selecting specific atoms participating in bonds
        cation_pi_selection = interactiontools.get_interaction_sel('cation_pi', group_a, group_b)
        pi_stacking_selection = interactiontools.get_interaction_sel('pi_stacking', group_a, group_b)
        #hbond_selection=interactiontools.get_interaction_sel('hbond', group_a,group_b)
        salt_bridge_selection = interactiontools.get_interaction_sel('salt_bridge', group_a, group_b)

    data = {
        'frame': [], 'n_cont': [], 'cont_a': [], 'cont_b': [], 'q_values': [] 
    }

    if bool_bonds:
        data['cation_pi'] = []
        data['pi_stacking'] = []
        data['hbond'] = []
        data['salt_bridge'] = []

    for ts in u.trajectory:
        frm = ts.frame
        dims = ts.dimensions
        frm_abs = k * chunk_size + frm

        cm_A = group_a.center_of_mass(wrap=True, unwrap=False, compound='group')
        cm_B = group_b.center_of_mass(wrap=True, unwrap=False, compound='group')
        dist_AB = contacts.distance_array(cm_A[None, :], cm_B[None, :], box=dims)[0, 0]
        
        if dist_AB < d_max * dims[2]:
            dist = contacts.distance_array(group_a, group_b, box=dims)
            mat_contacts = contacts.contact_matrix(dist, radius)
            n_contacts = mat_contacts.sum()

            if n_contacts > 0:
                pairs = np.array(np.where(mat_contacts))

                cont_a_ids_u = group_a.atoms.ids[pairs[0, :]]
                cont_b_ids_u = group_b.atoms.ids[pairs[1, :]]
                cont_a_rids_u = group_a.atoms.resids[pairs[0, :]]
                cont_b_rids_u = group_b.atoms.resids[pairs[1, :]]

                data['frame'].append(frm_abs)
                data['n_cont'].append(n_contacts)
                data['cont_a'].append(cont_a_ids_u.tolist())
                data['cont_b'].append(cont_b_ids_u.tolist())
                
                # Native contacts (q-values) calculation
                q_values_dict = {}
                if bool_q: 
                    pattern = os.path.join(output_path, "*_reference.pkl")
                    matching_files = glob.glob(pattern)
    
                    for file in matching_files:
                        ref_val = pd.read_pickle(file)
                        try:                        
                            ind_nat_A = np.asarray(ref_val["nc_a_ind"])
                            ind_nat_B = np.asarray(ref_val["nc_b_ind"])
                            r0 = np.asarray(ref_val["nc_dist"]).flatten()
                            prots = ref_val["prots"]
                            print("selections", len(group_a), len(group_b), flush=True)
        
                            if len(ind_nat_A) == 0 or len(ind_nat_B) == 0:
                                print(f"Skipping {prots}: empty reference contact set")
                                continue
                    
                            if np.max(ind_nat_A) >= dist.shape[0] or np.max(ind_nat_B) >= dist.shape[1]:
                                print(f"Skipping {prots}: indices out of bounds for current group")
                                continue
                                
                            r = dist[ind_nat_A, ind_nat_B]
                            q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))))
                            q_values_dict[prots] = float(q)
                        except:
                            print('no ' + ref_val["prots"] + ' interface')
                data['q_values'].append(q_values_dict)

                # Interaction profiles 
                if bool_bonds:
                    l1 = list(set(cont_a_rids_u))
                    l2 = list(set(cont_b_rids_u))
                    all_cont_resids = l1 + l2

                    all_chains = (group_a + group_b)
                    all_resids = all_chains.residues.resids

                    def eval_candidate_pairs(candiate_pairs, all_cont_resids, all_chains, all_resids):
                        """
                        Filter candidate atom pairs to ensure they belong to residues in contact.
                        
                        Parameters
                        ----------
                        candiate_pairs : array-like
                        all_cont_resids : list
                        all_chains : AtomGroup
                        all_resids : array-like
                        
                        Returns
                        -------
                        numpy.ndarray
                        """
                        mask_subset = np.isin(np.array(candiate_pairs), all_cont_resids).all(axis=1)
                        filtered_pairs = np.array(candiate_pairs)[mask_subset]

                        indices_a = np.searchsorted(all_resids, np.array(filtered_pairs)[:,0])
                        indices_b = np.searchsorted(all_resids, np.array(filtered_pairs)[:,1])
                        mols_a = all_chains.residues[indices_a]
                        mols_b = all_chains.residues[indices_b]
                        mol_pairs = np.column_stack((mols_a, mols_b))

                        return mol_pairs
                    
                    cation_pi_candidates = eval_candidate_pairs(cation_pi_selection, all_cont_resids, all_chains, all_resids)
                    cation_pi_contacts = interactiontools.cation_pi_contact(cation_pi_candidates, distance_cutoff=6.0, angle_cutoff=60.0)

                    pi_stacking_candidates = eval_candidate_pairs(pi_stacking_selection, all_cont_resids, all_chains, all_resids)
                    pi_stacking_contacts = interactiontools.pi_stacking_contact(pi_stacking_candidates, distance_cutoff=7.0, angle_cutoff=30.0, psi_cutoff=45.0)

                    salt_bridge_candidates = eval_candidate_pairs(salt_bridge_selection, all_cont_resids, all_chains, all_resids)
                    salt_bridge_contacts = interactiontools.salt_bridge_contact(salt_bridge_candidates, distance_cutoff=4.0)

               
                    #hbond #todo
                    #for donor -h
                    #mask_30 = np.isin(l_sel3[0].ids, cont_a_ids_u)

                    #for hydrogens -h

                    #atoms in contact (cont_a_ids_u) but from atom selection including hydrogens (u_group_a_sel)
                    #mask_30_c = np.isin(u_group_a_sel.ids, cont_a_ids_u)

                    #all atoms of the corresponding residues
                    #group_a_c=u_group_a_sel[mask_30_c].residues.atoms

                    #but only take the hydrogens from group_a_c
                    #group_c=group_a_c.select_atoms(group_a_sc_H)

                    #for aceptors -
                    #mask_31 = np.isin(l_sel3[1].ids, cont_b_ids_u)

                    
                    #hbond for later
                    #bool_hbond =interactiontools.hbond_contact(l_sel3[0][mask_30], group_c, l_sel3[1][mask_31], distance_cutoff=3.5, angle_cutoff=30.0)


                    data['cation_pi'].append(cation_pi_contacts)
                    data['pi_stacking'].append(pi_stacking_contacts)
                    data['hbond'].append([]) 
                    data['salt_bridge'].append(salt_bridge_contacts)
                    
    if bool_debug:
        print(f"[{worker_id}] RSS before trim: {psutil.Process().memory_info().rss / 1024 ** 2:.2f} MB", flush=True)

    #save memory
    del u, group_a, group_b, ts
    release_memory()
    
    if bool_debug:
        #print(f"{worker_id.memory_manager.memory_limit / 1024**3:.2f} GB", flush=True)
        print(f"[{worker_id}] memory after processing {traj}: {psutil.Process().memory_info().rss / 1024**2:.2f} MB", flush=True)

    
    if not data['frame']:
        return empty_contact_df()

    return pd.DataFrame(data)   


def contacts_within_cutoff_fast_parallel_split(top, group_a_sel, group_b_sel, traj_min, traj_lim, chunk_path, chunk_size, d_max=1.0, radius=4.5):
    """
    Generate delayed tasks for processing all trajectory chunks.
    
    Parameters
    ----------
    top : str
    group_a_sel, group_b_sel : str
    traj_min, traj_lim : int
    chunk_path : str
    chunk_size : int
    d_max : float
    radius : float
    
    Returns
    -------
    dask.dataframe.DataFrame
    """
    #l_chunk_files = sorted(glob.glob(chunk_path + "/chunk_*.xtc"))
    l_chunk_files = l_chunk_files = sorted(
        glob.glob(chunk_path + "/chunk_*.xtc"), 
        key=lambda x: int(re.search(r'chunk_(\d+)', x).group(1))
    )
    print('chunk_files',l_chunk_files, flush=True)

    delayed_chunks = [
        process_frame_chunk_for_split(top, traj, group_a_sel, group_b_sel, d_max, radius, k, chunk_size)
        for k, traj in enumerate(l_chunk_files)
    ]
    
    meta_dict = {
        'frame': pd.Series(dtype='int64'),
        'n_cont': pd.Series(dtype='int64'),
        'cont_a': pd.Series(dtype='object'),
        'cont_b': pd.Series(dtype='object'),
        'q_values': pd.Series(dtype='object'),
    }

    if bool_bonds:
        meta_dict.update({
            'cation_pi': pd.Series(dtype='object'),
            'pi_stacking': pd.Series(dtype='object'),
            'hbond': pd.Series(dtype='object'),
            'salt_bridge': pd.Series(dtype='object')
        })

    meta = pd.DataFrame(meta_dict)

    return dd.from_delayed(delayed_chunks, meta=meta)
        

@delayed
def calc_interactions(pair, top, traj, name, path, chunk_path, d_max, traj_min, traj_lim, delta_r, chunk_size, bool_sidechains=False, bool_bb=False):  
    """
    Orchestrate interaction calculation for a specific pair of molecules.
    
    Parameters
    ----------
    pair : tuple
    top : str
    traj : list
    name, path, chunk_path : str
    d_max, delta_r : float
    traj_min, traj_lim, chunk_size : int
    bool_sidechains, bool_bb : bool
    
    Returns
    -------
    Delayed object
    """
    job_id = os.environ.get('SLURM_JOB_ID', 'local_debug')
    tmp_path = f"/localscratch/{job_id}/tmp_parts"
    
    a = pair[0][0]
    b = pair[1][0]

    if bool_sidechains:
        A_mol_sel = f"({pair[0][1]}) and not name H* and not backbone and not name MW4"
        B_mol_sel = f"({pair[1][1]}) and not name H* and not backbone and not name MW4"

        A_mol_sel_h = f"({pair[0][1]}) and not backbone and not name MW4"
        B_mol_sel_h = f"({pair[1][1]}) and not backbone and not name MW4"

    elif bool_bb:
        A_mol_sel = f"({pair[0][1]}) and not name H* and backbone and not name MW4"
        B_mol_sel = f"({pair[1][1]}) and not name H* and backbone and not name MW4"

        A_mol_sel_h = f"({pair[0][1]}) and backbone and not name MW4"
        B_mol_sel_h = f"({pair[1][1]}) and backbone and not name MW4"
    
    else:
        A_mol_sel = f"({pair[0][1]}) and not name H* and not name MW4"
        B_mol_sel = f"({pair[1][1]}) and not name H* and not name MW4"

        A_mol_sel_h = f"({pair[0][1]}) and not name MW4"
        B_mol_sel_h = f"({pair[1][1]}) and not name MW4"

    print("a,b", a, b, flush=True)
    print("a_sel,b_sel", A_mol_sel_h, B_mol_sel_h, flush=True)

    contacts_ha_ddf = contacts_within_cutoff_fast_parallel_split(top, A_mol_sel_h, B_mol_sel_h, traj_min, traj_lim, chunk_path, chunk_size, d_max=d_max, radius=delta_r)
     
    return annotate_and_save_contacts(
        contacts_ha_ddf, A_mol_sel, B_mol_sel, top, traj, a, b, path, name)


def calc_contact_mat(top, traj, chains, dim, delta_r, path, chunk_size, chunk_path, name='cont_map', d_max=0.5, traj_min=0, traj_lim=5, bool_sidechains=False, bool_bb=False):
    """
    Generate the list of delayed tasks for all chain pairs.
    
    Parameters
    ----------
    top : str
    traj : list
    chains : list
    dim : array
    delta_r : float
    path : str
    chunk_size : int
    chunk_path : str
    name : str
    d_max : float
    traj_min, traj_lim : int
    bool_sidechains, bool_bb : bool
    
    Returns
    -------
    list of Delayed objects
    """
    delayed_path_tasks = [
        calc_interactions(pair, top, traj, name, path, chunk_path, d_max, traj_min, traj_lim, delta_r, chunk_size, bool_sidechains=bool_sidechains, bool_bb=bool_bb) 
        for pair in chains
    ]
    
    return delayed_path_tasks


# --- Helper Functions ---

def split_consecutive_intervals2(arr):
    """
    Splits a list of numbers into sublists of consecutive integers.
    
    Parameters
    ----------
    arr : list
    
    Returns
    -------
    list of lists
    """
    if not arr:
        return []

    arr = np.array(arr)
    split_indices = np.where(np.diff(arr) != 1)[0] + 1
    sublists = np.split(arr, split_indices)
    return [x.tolist() for x in sublists]


def get_matches_or_subsets(list_a, list_b):
    """
    Efficient O(N+M) matching of intervals.
    
    Parameters
    ----------
    list_a : list of lists
    list_b : list of lists
    
    Returns
    -------
    list
    """
    set_a_full = set(tuple(x) for x in list_a)
    set_b_lookup = set(tuple(x) for x in list_b)

    pool_a = [x for x in list_a if tuple(x) not in set_b_lookup]
    
    final_results = []
    pool_idx = 0
    n_pool = len(pool_a)

    for sub_b in list_b:
        tuple_b = tuple(sub_b)
        
        if tuple_b in set_a_full:
            final_results.append(sub_b)
            continue
            
        b_start = sub_b[0]
        b_end = sub_b[-1]
        
        while pool_idx < n_pool and pool_a[pool_idx][0] < b_start:
            pool_idx += 1
            
        temp_idx = pool_idx
        found_subset = None
        
        while temp_idx < n_pool:
            cand = pool_a[temp_idx]
            cand_start = cand[0]
            cand_end = cand[-1]
            
            if cand_start > b_end:
                break
            
            if cand_end <= b_end:
                found_subset = cand
                break 
            
            temp_idx += 1

        if found_subset:
            final_results.append(found_subset)

    return final_results


def process_overlap2(raw_list_a, raw_list_b):
    """
    Process overlaps between two raw lists of frame indices.
    
    Parameters
    ----------
    raw_list_a, raw_list_b : list
    
    Returns
    -------
    list
    """
    split_a = split_consecutive_intervals2(raw_list_a)
    split_b = split_consecutive_intervals2(raw_list_b)

    matches = get_matches_or_subsets(split_a, split_b)

    if not matches:
        return []
    
    return np.concatenate(matches).tolist()


def collect_list(arr):
    """
    Flatten list of sublists if sublist length is < eps_ts.
    
    Parameters
    ----------
    arr : list of lists
    
    Returns
    -------
    list
    """
    bla = [sublist for sublist in arr if len(sublist) < eps_ts]
    return [item for sublist in bla for item in sublist]


def split_consecutive_intervals(arr):
    """
    Return list of list according to intervals.
    
    Parameters
    ----------
    arr : list
    
    Returns
    -------
    list of lists
    """
    arr = np.array(arr)
    split_indices = np.where(np.diff(arr) != 1)[0] + 1
    sublists = np.split(arr, split_indices)
    
    if len(sublists) > 1:
        return [[sublists[0][0]], sublists[1].tolist()]
    else:
        return [sublists[0].tolist()]


def list_difference(list_a, list_b):
    """
    Return elements in list_a but not in list_b (set difference).
    
    Parameters
    ----------
    list_a, list_b : list
    
    Returns
    -------
    list
    """
    a = list_a if isinstance(list_a, list) else []
    b = list_b if isinstance(list_b, list) else []
    return list(np.setdiff1d(a, b))


def first_not_second(list_a, list_b):
    """
    Return elements in list_a but not in list_b (sorted list).
    
    Parameters
    ----------
    list_a, list_b : list
    
    Returns
    -------
    list
    """
    a = list_a if isinstance(list_a, list) else []
    b = list_b if isinstance(list_b, list) else []
    res = sorted(list(set(a) - set(b)))
    return res


def first_not_second_list(list_of_lists, list_b):
    """
    Return elements in each sublist of list_of_lists that are not in list_b, dropping empty sublists.
    
    Parameters
    ----------
    list_of_lists : list of lists
    list_b : list
    
    Returns
    -------
    list of lists
    """
    b = list_b if isinstance(list_b, list) else []
    if not isinstance(list_of_lists, list):
        return []
    
    return [sorted(list(set(sublist) - set(b))) for sublist in list_of_lists if set(sublist) - set(b)]


def first_and_second(list_a, list_b):
    """
    Return elements common to both list_a and list_b (intersection).
    
    Parameters
    ----------
    list_a, list_b : list
    
    Returns
    -------
    list
    """
    a = list_a if isinstance(list_a, list) else []
    b = list_b if isinstance(list_b, list) else []
    return sorted(list(set(a).intersection(set(b))))


def list_union(list_a, list_b):
    """
    Return the union of list_a and list_b.
    
    Parameters
    ----------
    list_a, list_b : list
    
    Returns
    -------
    list
    """
    a = list_a if isinstance(list_a, list) else []
    b = list_b if isinstance(list_b, list) else []
    return sorted(list(set(a).union(set(b))))
  
        
def load_contact_mats(parquet_path_ha, parquet_path_lb):
    """
    Load contact matrices from Parquet directories.
    
    Parameters
    ----------
    parquet_path_ha : str
    parquet_path_lb : str
    
    Returns
    -------
    tuple
    """
    contact_mat = pd.read_parquet(parquet_path_ha)
    try:
        contact_mat_lb = pd.read_parquet(parquet_path_lb)
    except:
        contact_mat_lb = None
    return contact_mat, contact_mat_lb


def atom_to_res(df):
    """
    Convert atom-level contacts to residue-level contacts.
    Aggregates information per frame and residue pair.
    
    Parameters
    ----------
    df : pandas.DataFrame
    
    Returns
    -------
    pandas.DataFrame
    """
    def frame_atom_lists(sub_df):
        """
        Aggregate contact info for a specific residue pair across frames.
        
        Parameters
        ----------
        sub_df : pandas.DataFrame
        
        Returns
        -------
        pandas.Series
        """
        frames = []
        q_values_list = []
        cation_pi_list = []
        pi_stacking_list = []
        hbond_list = []
        salt_bridge_list = []

        for frame, group in sub_df.groupby('frame'):
            frames.append(frame)
            
            if 'q_values' in group.columns and not group.empty:
                q_values_list.append(group['q_values'].iloc[0])
            else:
                q_values_list.append(None)
                
            if 'cation_pi' in group.columns and not group.empty:
                cation_pi_list.append(group['cation_pi'].iloc[0])
            else:
                cation_pi_list.append(None)

            if 'pi_stacking' in group.columns and not group.empty:
                pi_stacking_list.append(group['pi_stacking'].iloc[0])
            else:
                pi_stacking_list.append(None)
                
            if 'hbond' in group.columns and not group.empty:
                hbond_list.append(group['hbond'].iloc[0])
            else:
                hbond_list.append(None)

            if 'salt_bridge' in group.columns and not group.empty:
                salt_bridge_list.append(group['salt_bridge'].iloc[0])
            else:
                salt_bridge_list.append(None)
        
        return pd.Series({
            'frame': frames,
            'q_values': q_values_list,
            'cation_pi': cation_pi_list,
            'pi_stacking': pi_stacking_list,
            'hbond': hbond_list,
            'salt_bridge': salt_bridge_list
        })

    group_cols = ['res_a', 'res_b', 'type_a', 'type_b']
    result = (
        df.groupby(group_cols, group_keys=False)
          .apply(frame_atom_lists)
          .reset_index()
    )
    return result


def eval_contact_mat(res_contact_mat, res_contact_mat_lb, delta_ts=10, bool_breaks=True):
    """
    Evaluate contact matrix for valid contacts and breaks.
    
    Parameters
    ----------
    res_contact_mat : pandas.DataFrame
    res_contact_mat_lb : pandas.DataFrame
    delta_ts : int
    bool_breaks : bool
    
    Returns
    -------
    tuple
    """
    frames = res_contact_mat['frame']

    #check if there are contacts
    if len(frames) > 0:

        #identify breaks
        res_contact_mat['range_cont_tmp'] = list(res_contact_mat.apply(lambda row: list(np.arange(row['frame'][0], row['frame'][-1] + 1, 1)) if len(row['frame']) > 0 else [], axis=1))
        res_contact_mat['breaks'] = list(res_contact_mat.apply(lambda row: list_difference(row['range_cont_tmp'], row['frame']), axis=1))
        res_contact_mat.drop(columns=['range_cont_tmp'], inplace=True)
        df_breaks = res_contact_mat[res_contact_mat['breaks'].apply(lambda x: len(x) > 0)]
        
    else:

        #no breaks
        df_breaks = pd.DataFrame()
   
    if len(df_breaks) > 0 and bool_breaks:

        #analyse breaks - compare to lower bound cutoff (less strict)
        df_breaks_valid = pd.merge(
            df_breaks, 
            res_contact_mat_lb, 
            on=['res_a', 'res_b'], 
            how='inner', 
            suffixes=('_ha', '_lb')
        )

        df_breaks_valid = df_breaks_valid.drop(columns=['type_a_lb', 'type_b_lb'])
        df_breaks_valid = df_breaks_valid.rename(columns={
            'type_a_ha': 'type_a',
            'type_b_ha': 'type_b',
            'frame_ha': 'frame',
        })

        #candidates from lower bound cutoff (to be included)
        df_breaks_valid['frame_pot'] = list(df_breaks_valid.apply(lambda row: first_and_second(row['breaks'], row['frame_lb']), axis=1))

        #option a: identify and add frames from breaks: add frames from the break until the lower bound cutoff is exeeded 
        df_breaks_valid['frame_dp'] = [process_overlap2(pot, brk) for pot, brk in zip(df_breaks_valid['frame_pot'], df_breaks_valid['breaks'])]
        df_breaks_valid['frame_d2'] = list(df_breaks_valid.apply(lambda row: list_union(row['frame'], row['frame_dp']), axis=1))


        #option b: add all frames that are within the lower bound cutoff from the break 
        df_breaks_valid['frame_d'] = list(df_breaks_valid.apply(lambda row: list_union(row['frame'], first_and_second(row['breaks'], row['frame_lb'])), axis=1))

        df_breaks_valid['breaks_r'] = list(df_breaks_valid.apply(lambda row: first_not_second(row['breaks'], row['frame_d2']), axis=1).apply(split_consecutive_intervals))

        #join contacts (if eps_ts > 1, breaks are not longer than eps_ts)
        if eps_ts > 1:
            df_breaks_valid['frame_dtt'] = df_breaks_valid['breaks_r'].apply(collect_list)
            df_breaks_valid['frame_dt']  = list(df_breaks_valid.apply(lambda row: list_union(row['frame_d'], row['frame_dtt']), axis=1))
            df_breaks_valid['frame_d2t'] = list(df_breaks_valid.apply(lambda row: list_union(row['frame_d2'], row['frame_dtt']), axis=1))
        else:
            df_breaks_valid['frame_dt']  = df_breaks_valid['frame_d']
            df_breaks_valid['frame_d2t'] = df_breaks_valid['frame_d2']

        #remaining breaks
        #df_breaks_valid['breaks_rr'] = list(df_breaks_valid.apply(lambda row: first_not_second_list(row['breaks_r'], row['frame_dtt']), axis=1))
        res_contact_mat = pd.merge(res_contact_mat, df_breaks_valid[['res_a', 'res_b','frame_d','frame_lb', 'frame_dt','frame_d2t']], on=['res_a', 'res_b'], how='left')
        
        res_contact_mat[['frame_dt','frame_d','frame_lb','frame_d2t']] = res_contact_mat[['frame_dt','frame_d','frame_lb', 'frame_d2t']].fillna(0)
        
        # cleanup - deal with Nones and zeros
        clean_list = lambda x: [] if x == 0 else x
        res_contact_mat["frame_dt"]  = list(res_contact_mat.apply(lambda row: clean_list(row['frame_dt']), axis=1))
        res_contact_mat["frame_d2t"] = list(res_contact_mat.apply(lambda row: clean_list(row['frame_d2t']), axis=1))
        #res_contact_mat['breaks_rr'] = list(res_contact_mat.apply(lambda row: clean_list(row['breaks_rr']), axis=1))
        res_contact_mat['frame_d']   = list(res_contact_mat.apply(lambda row: clean_list(row['frame_d']), axis=1))
        res_contact_mat['frame_lb']  = list(res_contact_mat.apply(lambda row: clean_list(row['frame_lb']), axis=1))

        # use original values (in case no lower bound cutoffs are included)
        res_contact_mat["frame_dt"] = list(res_contact_mat.apply(lambda row: row["frame"] if len(row["frame"]) > len(row['frame_dt']) else row['frame_dt'], axis=1))
        res_contact_mat["frame_d2t"] = list(res_contact_mat.apply(lambda row: row["frame"] if len(row["frame"]) > len(row['frame_d2t']) else row['frame_d2t'], axis=1))
        
    elif len(frames) > 0 and not bool_breaks:
        print("contacts & no futher processing of breakes")
        res_contact_mat["frame_dt"]  = res_contact_mat["frame"]
        res_contact_mat["frame_d2t"] = res_contact_mat["frame"]
        res_contact_mat["frame_d"]   = np.nan
        res_contact_mat["frame_lb"]   = np.nan
        #res_contact_mat["breaks_rr"] = list(res_contact_mat.apply(lambda row: row['breaks'] if len(row['breaks']) > 0 else [], axis=1))
        return pd.DataFrame({'breaks': []}), res_contact_mat
    else: 
        print("No contacts & breaks to be detected")
        #res_contact_mat["breaks_rr"] = np.nan
        res_contact_mat["frame_dt"]  = np.nan
        res_contact_mat["frame_d2t"] = np.nan
        res_contact_mat["frame_d"]   = np.nan
        res_contact_mat["frame_lb"]   = np.nan
        return pd.DataFrame({'breaks': []}), res_contact_mat
    
    return df_breaks_valid, res_contact_mat


def pp_pair(pair, ha_result, lb_result, path, name, eps_ts, bool_breaks):
    """
    Post-process a single pair of molecules (worker function).
    
    Parameters
    ----------
    pair : tuple
    ha_result, lb_result : tuple
    path : str
    name : str
    eps_ts : int
    bool_breaks : bool
    
    Returns
    -------
    tuple
    """
    worker = get_worker()
    worker_id = worker.id if worker else "Unknown"

    a = pair[0][0]
    b = pair[1][0]
    print(f"[{worker_id}] Post-processing pair {a}-{b}", flush=True)

    lb_parquet_path = None
    contact_mat_lb = None

    _, ha_parquet_path = ha_result

    if ha_parquet_path is None:
        print(f"Skipping post-processing for pair {a}-{b} due to empty HA contact data.")
        return (a, b, "skipped")

    if bool_breaks:
        _, lb_parquet_path = lb_result
    
    contact_mat, contact_mat_lb = load_contact_mats(ha_parquet_path, lb_parquet_path)
    res_contact_matti = atom_to_res(contact_mat)

    if bool_breaks and lb_parquet_path is not None:
        res_contact_matti_lb = atom_to_res(contact_mat_lb)
        breaks_valid, dt_res_contact_mat = eval_contact_mat(res_contact_matti, res_contact_matti_lb, delta_ts=eps_ts, bool_breaks=bool_breaks)
        
        if bool_debug: 
	#possible intermediate results
            title_res = f'res_lb_{name}_{a}_{b}.parquet'
            res_contact_matti_lb.to_parquet(f'{path}/{title_res}', index=False)

            title_res_lb = f'r_lb_{name}_{a}_{b}.parquet'
            contact_mat_lb.to_parquet(f'{path}/{title_res_lb}', index=False)

    else:
        breaks_valid, dt_res_contact_mat = eval_contact_mat(res_contact_matti, None, delta_ts=eps_ts, bool_breaks=bool_breaks)

    
    #possible intermediate results
    
    if bool_debug:

        title_res = f'r_{name}_{a}_{b}.parquet'
        contact_mat.to_parquet(f'{path}/{title_res}', index=False)

        title_res = f'res_{name}_{a}_{b}.parquet'
        res_contact_matti.to_parquet(f'{path}/{title_res}', index=False)

    title_dt = f'dt_{name}_{a}_{b}.parquet'
    dt_res_contact_mat.to_parquet(f'{path}/{title_dt}', index=False)

    return (a, b, "processed")


def pp_contact_mat(chains, eps_ts, path, ha_tasks, lb_tasks, name='cont_map', bool_breaks=True):
    """
    Creates post-processing tasks by directly chaining them to their dependencies.
    
    Parameters
    ----------
    chains : list
    eps_ts : int
    path : str
    ha_tasks, lb_tasks : list
    name : str
    bool_breaks : bool
    
    Returns
    -------
    list
    """
    final_pp_tasks = []

    for i, pair in enumerate(chains):
        ha_result_task = ha_tasks[i]

        if bool_breaks:
            lb_result_task = lb_tasks[i]
        else:
            lb_result_task = None
         
        pp_task = dask.delayed(pp_pair)(
            pair, 
            ha_result_task,
            lb_result_task,
            path,  
            name, 
            eps_ts, 
            bool_breaks
        )
        final_pp_tasks.append(pp_task)
    
    return final_pp_tasks
    

def main(top, traj, bookkeeping_path, dim, delta_ha, combinations_in, eps_ts, time_status, traj_min, traj_max, bool_breaks, d_step, n_workers, split_parts, b_sidechains):
    """
    Main function to calculate and evaluate contact maps.
    
    Parameters
    ----------
    top : str
    traj : list
    bookkeeping_path : str
    dim : array
    delta_ha : float
    combinations_in : list
    eps_ts : int
    time_status : dict
    traj_min, traj_max : int
    bool_breaks : bool
    d_step, n_workers, split_parts : int
    b_sidechains : bool
    """
    # parallelization
    print("available resources (GB)", psutil.virtual_memory().total // 1024**3)
    
    start = traj_min
    cut = traj_max
    
    if bool_splid:
        chunk_folder_name = 'traj_parts_' + str(c0) + '_' + str(ck)
        chunk_folder = bookkeeping_path + "/" + chunk_folder_name
    else:
        chunk_folder = bookkeeping_path + "/traj_parts"

    chunk_size = split_trajectory_into_parts(top, traj, traj_min, traj_max, output_dir=chunk_folder, parts=split_parts, d_step=d_step)
    print("chunk size", chunk_size)
    
    n_pairs = len(combinations_in)
    print(n_pairs, "chain pairs to compute", flush=True)

    for pair in combinations_in:
        a = pair[0][0]
        b = pair[1][0]

        out_dir_ha = os.path.join(bookkeeping_path, f"cont_map_{a}_{b}")
        os.makedirs(out_dir_ha, exist_ok=True)

        out_dir_lb = os.path.join(bookkeeping_path, f"lb_cont_map_{a}_{b}")
        os.makedirs(out_dir_lb, exist_ok=True)
    print("Directory creation complete.")

    # Use .get() to avoid crashing if not running via SLURM
    job_id = os.environ.get('SLURM_JOB_ID', 'local_debug')
    local_dir = f"/localscratch/{job_id}"
    os.makedirs(local_dir, exist_ok=True)

    tmp_path = f"/localscratch/{job_id}/tmp_parts"
    os.makedirs(tmp_path, exist_ok=True)

    print("SLURM_JOB_ID:", job_id, local_dir)
    print("Running on node:", socket.gethostname(), flush=True)
    threads = 1

    cluster = LocalCluster(
        n_workers=n_workers, 
        threads_per_worker=threads,   
        local_directory=local_dir,
        memory_limit=mem_str,
        dashboard_address=":8788",
        processes=True
    )
    client = Client(cluster)

    print("Dask dashboard:", client.dashboard_link, flush=True) 
    print('hi paralellization', n_workers, split_parts, threads, flush=True)
    print(f"Hallo from Client: {client.scheduler_info()['id']}", flush=True)

    # Calculation for strict cutoff
    ha_tasks = calc_contact_mat(
        top, traj, combinations_in, dim, delta_ha, bookkeeping_path,
        chunk_size, chunk_folder, d_max=buffer, traj_min=start, traj_lim=cut,
        bool_sidechains=b_sidechains, bool_bb=b_bb
    )
    
    delta_l = delta_ha + delta_eps

    # Calculation for less strict cutoff
    if bool_breaks:
        lb_tasks = calc_contact_mat(
            top, traj, combinations_in, dim, delta_l, bookkeeping_path,
            chunk_size, chunk_folder, name='lb_cont_map', d_max=buffer, 
            traj_min=start, traj_lim=cut, bool_sidechains=b_sidechains, bool_bb=b_bb
        )
    else:
        lb_tasks = None

    # Post-processing
    pp_tasks = pp_contact_mat(combinations_in, eps_ts, bookkeeping_path, ha_tasks, lb_tasks, bool_breaks=bool_breaks)
    dask.visualize(pp_tasks, filename=bookkeeping_path + '/full_workflow_graph.svg')

    # Dask report
    report = performance_report(filename=bookkeeping_path + "/dask-report.html")
    with report:
        pp_results = client.compute(pp_tasks, sync=True)
    
    print("Post-processing finished.", flush=True)

    if bool_splid:
        cleanup_temp_folders(bookkeeping_path, name=chunk_folder_name, complete=False)
    else:
        cleanup_temp_folders(bookkeeping_path, complete=True)


if __name__ == "__main__":   
    
    bookkeeping_path = output_path
    os.makedirs(bookkeeping_path, exist_ok=True)
    
    # Load domains
    domains = pd.read_parquet(f'{path}/{label}_domains.parquet')
    sys_domains = pd.read_parquet(f'{path}/{label}_sys_domains.parquet')

    # Load Trajectory data
    traj = [f'{path}/{label}_full_pi_ref.xtc']
    top1 = f"{path}/{label}_pi_clean.pdb"
    top2 = f"{path}/{label}_pi_clean.pdb"

    # If different data for ion analysis is used
    if bool_pwi:
        traj = [f'{path}/{label}_full_pi_ref.xtc']
        top1 = f"{path}/{label}_pi_clean.pdb"
        top2 = f"{path}/{label}_pi_clean.pdb"
        print(traj, len(traj), flush=True)
    
    # Check which file exists
    if os.path.exists(top2):
        top = top2
        u = mda.Universe(top, traj, continuous=True, in_memory=False)
    elif os.path.exists(top1):
        top = top1
        u = mda.Universe(top, traj, continuous=True, in_memory=False)
    else:
        raise FileNotFoundError("Neither topology file found.")
    print("using", top, "as topology")


    # Select ions and molecules
    ALL_MOLS_SEL = [f"resid {row['min']}-{row['max']}" for _, row in sys_domains.iterrows()]
    I_SEL = ["resname CL or resname NA"] #ions and water

    #check water amount and form batches
    w_call="resname SOL"
    #w_all=u.select_atoms(w_call)
    #w_batch_size = 2000
    #w_residues = u.select_atoms(w_call).residues #u1.residues[u1.residues.resnames=="SOL"]
    # make list of residue groups
    #res_batches = [w_residues[i:i + w_batch_size] for i in range(0, len(w_residues), w_batch_size)]
    #sol_batches = [batch.atoms for batch in res_batches]
    #W_SEL = [f"resid {' '.join(str(r.resid) for r in batch)}" for batch in res_batches]
    #print("water devided into ",len(W_SEL), " batches of ",w_batch_size)

    if bool_splid: 
        combinations_in = list(itertools.combinations(enumerate(ALL_MOLS_SEL), 2))[c0:ck] 
    elif bool_debug:
        combinations_in = list(itertools.combinations(enumerate(ALL_MOLS_SEL), 2))[10:11]
    else: 
        combinations_in = list(itertools.combinations(enumerate(ALL_MOLS_SEL), 2)) 
    

    # mod: skip chain pairs that have already been calculated in a previous run
    if bool_filtering:
        print("original combinations", combinations_in[0:10])
        pattern_done = "dt_cont_map_*_*.parquet"
        md_files = glob.glob(os.path.join(bookkeeping_path, pattern_done))
        found_ids = []
        regex = re.compile(r"dt_cont_map_(\d+)_(\d+)\.parquet")
        
        for file_path in md_files:
            filename = os.path.basename(file_path) 
            match = regex.search(filename)
            if match:
                id1 = int(match.group(1))
                id2 = int(match.group(2))
                found_ids.append((id1, id2))
                
        print(f"Already done: {len(found_ids)}")
        exclusion_set = set(found_ids)
        filtered_data = []
    
        for entry in combinations_in:      
            index_pair = (entry[0][0], entry[1][0])
            if index_pair not in exclusion_set:
                filtered_data.append(entry)
        
        print(f"Original length: {len(combinations_in)}")
        print(f"Filtered length: {len(filtered_data)}")
        print("filtered combinations", filtered_data[0:10])
    
    print(len(combinations_in), "from possible combinations", len(list(itertools.combinations(enumerate(ALL_MOLS_SEL), 2))), flush=True)

    # run
    combinations_in_ions = [
        ((i, a), (len(ALL_MOLS_SEL) + j, b))
        for i, a in enumerate(ALL_MOLS_SEL)
        for j, b in enumerate(I_SEL)
    ]

    if bool_pwi & bool_filtering:
        combinations = filtered_data + combinations_in_ions
        print("take filtered ion data")
    elif bool_filtering:
        print("take filtered data")
        combinations = filtered_data
    elif bool_pwi:
        combinations = combinations_in + combinations_in_ions
        print("take ion data")
    else:
        combinations = combinations_in
        print("take mol data")

    if traj_max == 0:
        traj_max = len(u.trajectory) - 1

    t_min = u.trajectory[traj_min].time
    t_max = u.trajectory[traj_max].time

    if bool_debug:
        traj_max = 2000
    
    # trajectory parameters
    ts = (u.trajectory[-1].time - u.trajectory[-2].time) * d_step 
    frames = int(len(u.trajectory[traj_min:traj_max]) / d_step)
    time_status = dict({'t_0': t_min, 't_max': t_max, 'ts': ts, 'nframes': frames, 'eps_ts': eps_ts})
    pd.DataFrame([time_status]).to_parquet(bookkeeping_path + '/' + label + '_time_status.parquet', index=False)
    print("contacts from", t_min, "to", t_max, "in", ts, "ps steps")

    contact_status = dict({'buffer': buffer, 'delta_ha': delta_ha, 'delta_eps': delta_eps, 'bool_breaks': bool_breaks})
    pd.DataFrame([contact_status]).to_parquet(bookkeeping_path + '/' + label + '_contact_status.parquet', index=False)

    main(top, traj, bookkeeping_path, u.dimensions, delta_ha, combinations, eps_ts, time_status, traj_min, traj_max, bool_breaks, d_step, n_workers, split_parts, b_sidechains)
