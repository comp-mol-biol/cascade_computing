import pandas as pd
import numpy as np
import glob
import sys
import os
import MDAnalysis as mda
import argparse
import pyarrow
from distutils.util import strtobool
import ast
import re
import json
from pandarallel import pandarallel


NP_INT_PATTERN = re.compile(r'np\.int64\((\d+)\)')

def get_args(argv=None):

    parser = argparse.ArgumentParser(description='input_output')
    parser.add_argument('--path_tools', action="store", dest='path_tools', default='./')
    parser.add_argument('--path', action="store", dest='path', default='./')
    parser.add_argument('--output', action="store", dest='output', default='./')
    parser.add_argument('--label', action="store", dest='label', default='ohoh')
    parser.add_argument('--result', action="store", dest='results_name', default='result_files')
    parser.add_argument('--pwi', action="store", dest='bool_pwi', default='False')
    parser.add_argument('--bonds', action="store", dest='bool_bonds', default='False') 
    parser.add_argument('--n_workers', action="store", dest='n_workers', default=64) 

    return parser.parse_args(argv)


def clean_and_filter_combined(row, key):
    """
    Parses, cleans, and filters in a single pass.
    """
    data = row[key]
    
    # 1. Access pandas/row values once (outside the loop)
    res_a = row['res_a']
    res_b = row['res_b']

    # 2. Handle top-level container types
    if isinstance(data, np.ndarray):
        data = data.tolist()
    if not isinstance(data, list):
        return []

    cleaned_result = []

    # 3. Iterate once through the outer list
    for sublist in data:
        # A. Clean the sublist if it is a string
        if isinstance(sublist, str):
            # Regex replace np.int64(...) -> ...
            processed_str = NP_INT_PATTERN.sub(r'\1', sublist)
            try:
                # Fast JSON parse (handling single quotes)
                sublist = json.loads(processed_str.replace("'", '"'))
            except (json.JSONDecodeError, TypeError):
                # If parsing fails, we can't iterate over it, so treat as empty
                cleaned_result.append([])
                continue

        # B. Filter Logic (Nested Loop)
        matches = []
        if isinstance(sublist, list):
            for item in sublist:
                # Recursive cleaning: If item is a string, parse it too
                if isinstance(item, str):
                    processed_item = NP_INT_PATTERN.sub(r'\1', item)
                    try:
                        item = json.loads(processed_item.replace("'", '"'))
                    except (json.JSONDecodeError, TypeError):
                        pass # Keep original item if parse fails

                # Ensure item is a list before checking contents
                if isinstance(item, list):
                    # C. The logic from 'parse_and_clean': convert digits to ints
                    # We do this JIT (Just-In-Time) before comparison
                    item = [
                        int(x) if (isinstance(x, str) and x.isdigit()) else x 
                        for x in item
                    ]

                    # D. The logic from 'filter_by_res': Check conditions
                    if len(item) == 2:
                        # Compare against cached row values
                        if item[0] == res_a and item[1] == res_b:
                            matches.append(item)
                
                # Note: We silently ignore non-list items (rogue integers) here,
                # effectively filtering them out.

        cleaned_result.append(matches)

    return cleaned_result



def split_consecutive_intervals(arr):
    """Split array into sublists of consecutive integers, ignoring NaNs.
    Always return at least one (possibly empty) list."""
    # Convert to numpy array
    arr = np.array(arr, dtype=float)  # Allows for NaNs

    # Remove NaNs
    arr = arr[~np.isnan(arr)]

    # If empty after cleaning, return one empty list
    if arr.size == 0:
        return [[]]

    # Convert to int
    arr = arr.astype(int)

    # If only one element, return it as a sublist
    if arr.size == 1:
        return [[int(arr[0])]]

    # Find breaks in consecutive integers
    split_indices = np.where(np.diff(arr) != 1)[0] + 1
    sublists = np.split(arr, split_indices)

    # Convert to plain Python ints
    return [[int(x) for x in sublist] for sublist in sublists]

def process_all_interactions(row):
        return pd.Series({
            'cation_pi_filtered': clean_and_filter_combined(row, 'cation_pi'),
            'pi_stacking_filtered': clean_and_filter_combined(row, 'pi_stacking'),
            'hbond_filtered': clean_and_filter_combined(row, 'hbond'),
            'salt_bridge_filtered': clean_and_filter_combined(row, 'salt_bridge')
        })


def symmetrize_by_id_pair(asymmetric_df, id_col_a, id_col_b):
        """
        Symmetrizes a DataFrame based on a primary ID pair, automatically swapping all related columns.

        Args:
            asymmetric_df (pd.DataFrame): DataFrame with one-way contacts.
            id_col_a (str): The main 'A' side identifier (e.g., 'res_a').
            id_col_b (str): The main 'B' side identifier (e.g., 'res_b').

        Returns:
            pd.DataFrame: A new DataFrame with fully symmetric contact pairs.
        """
        # 1. Find all corresponding '_a' and '_b' columns automatically
        cols_a = sorted([col for col in asymmetric_df.columns if col.endswith('_a')])
        cols_b = sorted([col for col in asymmetric_df.columns if col.endswith('_b')])

        if id_col_a not in cols_a or id_col_b not in cols_b:
            raise ValueError(f"'{id_col_a}' or '{id_col_b}' not found in DataFrame columns with _a/_b suffixes.")

        if len(cols_a) != len(cols_b):
            raise ValueError("Mismatched number of '_a' and '_b' columns found.")

        # 2. Create the rename map for all found pairs
        swap_map = dict(zip(cols_a + cols_b, cols_b + cols_a))

        # 3. Create the swapped copy
        df_swapped = asymmetric_df.rename(columns=swap_map)

        # 4. Combine and return
        symmetric_df = pd.concat([asymmetric_df, df_swapped], ignore_index=True)

        return symmetric_df

def main(path,output_path,contact_path,label,results_name,bool_pwi,bool_bonds):

    

    #################### load system ##########################
    top1 = f"{path}/{label}_pi_clean.pdb"

    if bool_pwi:
        top1=f"{path}/{label}_pi_clean.pdb"

    #load universe
    u = mda.Universe(top1)

    #time data
    time_df = pd.read_parquet(contact_path+'/'+label+'_time_status.parquet')
    n_frames=time_df['nframes'][0] 
    
    #load contact metadata
    pattern = os.path.join(path, "*_job_dict.parquet")
    matching_files = glob.glob(pattern)
    
    if matching_files:
        job_dict = pd.read_parquet(matching_files[0])

    try:
        all_domain_names=job_dict['domains_ordered'].iloc[0].split(";")
    except:
        all_domain_names=[]
        
    # load domains
    domains = pd.read_parquet(f'{path}/{label}_domains.parquet')
    sys_domains= pd.read_parquet(f'{path}/{label}_sys_domains.parquet')

    data_path=contact_path+"/"+results_name
    #data_path=job_dict['results_path'].iloc[0] #contact_path+'/'+results_name
    print("data_path", data_path, flush=True)
    os.makedirs(data_path, exist_ok=True)
    

    #################### evaluate meta data ##########################
    sys_domains['mol_ind'] = sys_domains.index
    try:
        sys_domains["struc_id"]=sys_domains['ind_mol']
    except:
        pass

    if "struc_id" not in sys_domains.columns:
        sys_domains["struc_id"]=np.nan


    # add original residue indices
    for i, row in sys_domains.iterrows():
        match = domains[domains["prot"] == row["prot"]]

        if not match.empty:  # Ensure match exists before accessing values
            min_val = match["min"].values[0]
            max_val = match["max"].values[0]

            print(i,[min_val, max_val])  # Debugging output

            # Assign list to the column using .at[]
            sys_domains.at[i, "org_min"] = min_val
            sys_domains.at[i, "org_max"] = max_val
    
   
    #################### meta data lookup table ##########################

    domain_names = [
    dom for dom in all_domain_names
    if f"{dom}_min" in sys_domains.columns and f"{dom}_max" in sys_domains.columns
]
    domain_intervals = {
    dom: pd.IntervalIndex.from_arrays(sys_domains[f"{dom}_min"], sys_domains[f"{dom}_max"], closed='both')
    for dom in domain_names
    if f"{dom}_min" in sys_domains and f"{dom}_max" in sys_domains
}

    residues = np.arange(sys_domains['min'].min(), sys_domains['max'].max() + 1)
    mol_intervals = pd.IntervalIndex.from_arrays(sys_domains['min'], sys_domains['max'], closed='both')

     
    res_df = pd.DataFrame({'res': residues})
    
    # lookup intervals for molecules and domains
    res_df['mol_ind'] = res_df['res'].apply(lambda x: mol_intervals.get_indexer([x])[0])
    res_df['res_dom']= np.nan 
    for dom in domain_names:
        res_df['res_dom'] = res_df.apply(
        lambda row: dom if domain_intervals[dom].get_indexer([row['res']])[0] > -1 else row['res_dom'], 
        axis=1)
 
    # Retrieve residue names from `uni`
    res_df['res_type'] = res_df['res'].apply(lambda x: u.residues[x - 1].resname)

    # add moltype from domains_sys
    res_df = res_df.merge(
        sys_domains[['mol_ind','struc_id','prot']], left_on='mol_ind', right_on='mol_ind', how='left'
    )

    #add original indicies
    n_mols=res_df['mol_ind'].values.max()
    for n in np.arange(n_mols+1):
        match=res_df[res_df['mol_ind']==n].copy() 
        ind_min=match.index.min()
        ind_max=match.index.max()
        #print(ind_min,ind_max)
        ind_org=match['res']-ind_min
        #print(len(ind_org))
        match = match.merge(
        domains[['prot','min','max']], left_on='prot', right_on='prot', how='left')
        val=ind_org.values+match['min']-1
        res_df.loc[match['res']-1, 'res_org'] = val.values
    res_df['res_org']=res_df['res_org'].fillna(-1).astype(int)

    if bool_pwi:
        #add ions/solvent
    
        ions_rg = u.select_atoms("resname CL or resname NA").residues
        ions_df = pd.DataFrame({ "res": ions_rg.resids.astype(int),
            "prot": str("ION"),
            "mol_ind": None,
            "res_org": None,
            "res_type": ions_rg.resnames})
        res_df = pd.concat([res_df, ions_df], ignore_index=True)

    # save RAM
    int_cols = ['res', 'mol_ind', 'struc_id', 'res_org']
    for col in int_cols:
        # Downcast to the smallest possible integer type
        res_df[col] = pd.to_numeric(res_df[col], downcast='integer')

    # meta data types
    category_cols = ['res_dom', 'res_type', 'prot']
    for col in category_cols:
        res_df[col] = res_df[col].astype('category')

    # save meta dataframe as lookup table
    title = label+'_contacts_res_meta'
    res_df.to_parquet(contact_path + '/' + title + '.parquet',engine="pyarrow")

    print("contact_path",contact_path)

    ####################load & cleanup contact data ##########################
    
    pattern = '*dt_cont_map*.parquet'
    files = glob.glob(os.path.join(contact_path, pattern))

    # Load and combine contact map data from all files
    data_comb = []

    bool_pp=True
    if bool_pp==False:
        for f in files:
            df = pd.read_parquet(f,columns=['res_a','res_b','frame', 'frame_dt','frame_lb'])
            data_comb.append(df)
    else: 
        for f in files:
            #df = pd.read_parquet(f,columns=['res_a','res_b','frame', 'frame_dt','breaks_rr', 'added_lb_frames', 'added_gap_frames','cation_pi','pi_stacking','hbond','salt_bridge','cation_pi_dt','pi_stacking_dt','hbond_dt','salt_bridge_dt'])
            try:
                df = pd.read_parquet(f,columns=['res_a','res_b','frame','frame_lb', 'frame_dt','frame_d2t','cation_pi','pi_stacking','hbond','salt_bridge'])
            except:
                df = pd.read_parquet(f,columns=['res_a','res_b','frame','frame_lb', 'frame_dt','frame_d2t'])
                print('cool no bonds')
                pass
            #make more data efficient
            data_comb.append(df)

    # concatenate all contact map data into a single DataFrame
    if len(data_comb)>1:
        df_contact_map_all = pd.concat(data_comb, keys=np.arange(len(data_comb)))

    else:
        df_contact_map_all =data_comb[0]
    print("df_contact_map_all, loaded", flush=True)

    #add length of frames in contact
    #df_contact_map_all['len_contacts'] = df_contact_map_all['frame_dt'].apply(
    #lambda x: len(x) if isinstance(x, (list, tuple, str, np.ndarray)) else 0)

    df_contact_map_all['len_contacts'] = df_contact_map_all['frame_d2t'].apply(
    lambda x: len(x) if isinstance(x, (list, tuple, str, np.ndarray)) else 0)


    # save RAM
    contact_int_cols = ['res_a', 'res_b'] # Add other integer columns if they exist
    for col in contact_int_cols:
        df_contact_map_all[col] = pd.to_numeric(df_contact_map_all[col], downcast='integer')

    print("Optimized contact map memory usage:", df_contact_map_all.info(memory_usage='deep'))


    #add meta data
    df_contact_map_all = df_contact_map_all.merge(
    res_df,
    left_on='res_a',
    right_on='res',
    how='left',
).drop(columns=['res'], errors='ignore')

    df_contact_map_all = df_contact_map_all.merge(
    res_df,
    left_on='res_b',
    right_on='res',
    how='left',
    suffixes=('_a', '_b')  # suffix for columns from res_df
).drop(columns=['res'], errors='ignore')

    print("pre-symmetry done", flush=True)
    #add symmetry - a/b values also added as b/a w.r.t. residue IDs 
    dt_cont_map_all_sym=symmetrize_by_id_pair(df_contact_map_all, 'res_a', 'res_b')
    print("symmetry done", flush=True)
    #save to file 
    title='_contacts_only'
    
    
    ####################percistance times############################
    #################################################################
   
    dt_cont_map_all_sym["percistance_times"]=dt_cont_map_all_sym["frame_d2t"].parallel_apply(lambda x: [len(sub) for sub in split_consecutive_intervals(x)])
    dt_cont_map_all_sym["avg_persistence"] = dt_cont_map_all_sym["percistance_times"].parallel_apply(lambda x: np.mean(x) if x else 0 )
    dt_cont_map_all_sym["median_persistence"] = dt_cont_map_all_sym["percistance_times"].parallel_apply( lambda x: np.median(x) if x else 0 )

    print(bool_bonds, '----------------------------------------------------------------------------------')
    if bool_bonds==True:

        ## process_all_interactions
        new_columns = dt_cont_map_all_sym.parallel_apply(process_all_interactions, axis=1)
        dt_cont_map_all_sym = pd.concat([dt_cont_map_all_sym, new_columns], axis=1)

        print("filtering Done.")
        
        ## include frame lists of lists
        frames_list = dt_cont_map_all_sym['frame'].tolist()
        
        cols_to_process = [
            ('cation_pi_filtered', 'cation_pi_t'),
            ('pi_stacking_filtered', 'pi_stacking_t'),
            ('hbond_filtered', 'hbond_t'),
            ('salt_bridge_filtered', 'salt_bridge_t')
        ]
        
        for source_col, target_col in cols_to_process:
            # Convert data column to list of lists
            data_list = dt_cont_map_all_sym[source_col].tolist()
            
            # Run the fast list comprehension (Zipping two lists together)
            dt_cont_map_all_sym[target_col] = [
                # Note: I simplified your logic to "Keep 'a' if 'b' is not empty"
                # because 'a[len(b)>0][0]' is often unstable/error-prone.
                [a for a, b in zip(row_frames, row_data) if len(b) > 0]
                for row_frames, row_data in zip(frames_list, data_list)
            ]
        
        print("pp Done.")

    else:
        dt_cont_map_all_sym['cation_pi_t']=dt_cont_map_all_sym['cation_pi']
        dt_cont_map_all_sym['pi_stacking_t']=dt_cont_map_all_sym['pi_stacking']
        dt_cont_map_all_sym['hbond_t']=dt_cont_map_all_sym['hbond']
        dt_cont_map_all_sym['salt_bridge_t']=dt_cont_map_all_sym['salt_bridge']
        print('no_cleaning_for_bond_info')

    try:
        if bool_bonds==True:
            print('saving with special bonds')
            clmns=['res_a', 'res_b', 'frame','frame_lb','frame_dt','frame_d2t','len_contacts','mol_ind_a', 'res_dom_a', 'res_type_a', 'struc_id_a', 'prot_a','res_org_a', 'mol_ind_b', 'res_dom_b', 'res_type_b', 'struc_id_b','prot_b', 'res_org_b','percistance_times','avg_persistence', 'median_persistence', 'cation_pi_t', 'pi_stacking_t', 'hbond_t', 'salt_bridge_t']
            dt_cont_map_all_sym[clmns].to_parquet(
            f'{data_path}/{label}{title}_t.parquet',
            engine='pyarrow',              # best support for column-wise reads
            compression='snappy',          # good balance of speed & compression
            index=False                    # skip index to save space
        )
        else:
            print('saving without special bonds')
            clmns=['res_a', 'res_b', 'frame','frame_lb','frame_dt','frame_d2t','len_contacts','mol_ind_a', 'res_dom_a', 'res_type_a', 'struc_id_a', 'prot_a','res_org_a', 'mol_ind_b', 'res_dom_b', 'res_type_b', 'struc_id_b','prot_b', 'res_org_b','percistance_times','avg_persistence', 'median_persistence']
            dt_cont_map_all_sym[clmns].to_parquet(
            f'{data_path}/{label}{title}_t.parquet',
            engine='pyarrow',              # best support for column-wise reads
            compression='snappy',          # good balance of speed & compression
            index=False                    # skip index to save space
        )
    except:
        print("saving exception", flush=True)
        dt_cont_map_all_sym[['res_a', 'res_b', 'frame_dt','len_contacts',
       'mol_ind_a', 'res_dom_a', 'res_type_a', 'struc_id_a', 'prot_a',
       'res_org_a', 'mol_ind_b', 'res_dom_b', 'res_type_b', 'struc_id_b',
       'prot_b', 'res_org_b','percistance_times','avg_persistence', 'median_persistence']].to_parquet(
        f'{data_path}/{label}{title}_t.parquet',
    engine='pyarrow',              # best support for column-wise reads
    compression='snappy',          # good balance of speed & compression
    index=False                    # skip index to save space
)
    
    return 0



if __name__ == "__main__":

    #parser
    args = get_args()

    n_workers = int(args.n_workers)

    #initialize parallel CPUs
    pandarallel.initialize(nb_workers=n_workers, use_memory_fs=True,progress_bar=True)
    

    main( path = args.path,output_path = args.output,contact_path=args.output,label=args.label,results_name= args.results_name,
    bool_pwi= bool(strtobool(args.bool_pwi)),
    bool_bonds=bool(strtobool(args.bool_bonds)))
