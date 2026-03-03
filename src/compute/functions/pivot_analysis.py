import pandas as pd
import numpy as np
import glob
import sys
import os
import pickle
import MDAnalysis as mda
from tqdm import tqdm
import argparse
from itertools import chain
from collections import Counter
import pyarrow
import itertools
from distutils.util import strtobool
import seaborn as sbn
import ast
import re
import json

from ..utils import toolbox_pp as pptools

amino_acids = [
'ALA', 'ARG', 'ASN', 'ASP', 'CYS',
'GLN','GLU', 'GLY', 'HIS', 'ILE',
'LEU', 'LYS', 'MET', 'PHE', 'PRO',
'SER', 'THR', 'TRP', 'TYR', 'VAL']


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

    return parser.parse_args(argv)


#evaluation functions
def union_length(series):
    return len(set(chain.from_iterable(series.dropna())))
                                       
def union(series):
    return list(set(chain.from_iterable(series.dropna())))

def make_list(series):
    return list(chain.from_iterable(series.dropna()))

def make_listol(series):
    return series.dropna().tolist()


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


def analyse_attribute(attr, attr_val, dt_cont_map_all_sym, data_path,args):

        path = args.path
        path_tools = args.path_tools
        output_path = args.output
        contact_path=output_path
        label=args.label
        results_name= args.results_name
        bool_pwi= bool(strtobool(args.bool_pwi))
        bool_bonds=bool(strtobool(args.bool_bonds))
    

        #################### available pairs (in system) ##########################
        
        #residue pairs - 
        grouped = dt_cont_map_all_sym.groupby([f'{attr}_a', f'{attr}_b'],sort=False)
        unique_pairs = grouped.apply(
            lambda df: len(df[['res_a', 'res_b']].drop_duplicates().to_records(index=False).tolist()))
        
        attr_unique_res_count=unique_pairs.unstack(fill_value=0)
        attr_unique_res_count = attr_unique_res_count.reindex(index=attr_val, columns=attr_val, fill_value=0)
        

        #orig residue pairs
        unique_pairs = grouped.apply(
            lambda df: len(df[['res_org_a', 'res_org_b']].drop_duplicates().to_records(index=False).tolist()))
        
        attr_unique_res_org_count=unique_pairs.unstack(fill_value=0)
        attr_unique_res_org_count = attr_unique_res_org_count.reindex(index=attr_val, columns=attr_val, fill_value=0)


        #chain pairs
        unique_pairs = grouped.apply(
            lambda df: len(df[['mol_ind_a', 'mol_ind_b']].apply(lambda row: tuple(sorted((row['mol_ind_a'], row['mol_ind_b']))), axis=1).drop_duplicates().tolist()))

        attr_unique_chains_count=unique_pairs.unstack(fill_value=0)
        attr_unique_chains_count = attr_unique_chains_count.reindex(index=attr_val, columns=attr_val, fill_value=0)


        #residue type pairs
        unique_pairs = grouped.apply(
            lambda df: len(df[['res_type_a', 'res_type_b']].drop_duplicates().to_records(index=False).tolist()))
        
        attr_unique_res_type_count=unique_pairs.unstack(fill_value=0)
        attr_unique_res_type_count = attr_unique_res_type_count.reindex(index=attr_val, columns=attr_val, fill_value=0)


        #save to file
        attr_unique_res_count.to_parquet(f"{data_path}/{label}_{attr}_count_res.parquet",engine="pyarrow")
        attr_unique_res_org_count.to_parquet(f"{data_path}/{label}_{attr}_count_res_orig.parquet",engine="pyarrow")
        attr_unique_chains_count.to_parquet(f"{data_path}/{label}_{attr}_count_chains.parquet",engine="pyarrow")
        attr_unique_res_type_count.to_parquet(f"{data_path}/{label}_{attr}_count_res_type.parquet",engine="pyarrow")


        #################### pivot tables ##########################

        
        # len(C(i,j)) aggregated (union of all unique contact frames)
        contacts_pivot_prob = dt_cont_map_all_sym.pivot_table(
        index=attr+"_a",
        columns=attr+"_b",
        values='frame_d2t',
        aggfunc=union_length,
        observed=False
    )
        contacts_pivot_full_prob = contacts_pivot_prob.reindex(index=attr_val, columns=attr_val,fill_value=0) 
        title='_'+attr+'_prob_pivot.parquet'
        contacts_pivot_full_prob.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        #contacts_pivot_full_prob.to_pickle(f'{data_path}/{label}_{attr}_prob_pivot.pkl')


        # L(i,j) (sum of all individual contact frames)
        contacts_pivot_time = dt_cont_map_all_sym.pivot_table(
        index=attr+"_a",
        columns=attr+"_b",
        values='len_contacts',
        aggfunc='sum',
        observed=False
    )
        contacts_pivot_full_time = contacts_pivot_time.reindex(index=attr_val, columns=attr_val,fill_value=0)
        title='_'+attr+'_time_count_pivot.parquet'
        contacts_pivot_full_time.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        #contacts_pivot_full_time.to_pickle(f'{data_path}/{label}_{attr}_time_count_pivot.pkl')


        # C(i,j) aggregated (union of all unique contact frames)
        distr_pivot_frame_union = dt_cont_map_all_sym.pivot_table(
        index=attr+"_a",
        columns=attr+"_b",
        values='frame_d2t',
        aggfunc=union,
        observed=False
        )

        distr_pivot_frame_union = distr_pivot_frame_union.reindex(index=attr_val, columns=attr_val,fill_value=0)
        distr_pivot_frame_union=distr_pivot_frame_union.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))

        title='_'+attr+'_frame_union_pivot.parquet'
        distr_pivot_frame_union.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        #title='_'+attr+'_frame2_union_pivot.pkl'
        #distr_pivot_frame_union2.to_pickle(f'{data_path}/{label}{title}')
        

        # C(i,j) aggregated ---> distribution of percistance
        #distr_pivot_percistance=distr_pivot_frame_union.applymap(
        #lambda x: [len(sub) for sub in split_consecutive_intervals(x)] if isinstance(x, list) else x).fillna(0)
        distr_pivot_percistance = distr_pivot_frame_union.applymap(lambda x: [len(sub) for sub in split_consecutive_intervals(x)]
    if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else [])
)

        distr_pivot_percistance_full = distr_pivot_percistance.reindex(index=attr_val, columns=attr_val,fill_value=0)
        distr_pivot_percistance_full=distr_pivot_percistance_full.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))
                          
        title='_'+attr+'_agg_percistance_distribution_pivot.parquet'
        distr_pivot_percistance_full.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        #title='_'+attr+'_agg_percistance_distribution_pivot.pkl'
        #distr_pivot_percistance_full.to_pickle(f'{data_path}/{label}{title}')

        
        #list[C_k(i,j)] (list of all unique contact frames)
        contacts_pivot_list = dt_cont_map_all_sym.pivot_table(
        index=attr+"_a",
        columns=attr+"_b",
        values='frame_d2t',
        aggfunc=make_listol,
        observed=False
    )
        contacts_pivot_list=contacts_pivot_list.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))
        contacts_pivot_full_list = contacts_pivot_list.reindex(index=attr_val, columns=attr_val,fill_value=0) 
        contacts_pivot_full_list=contacts_pivot_full_list.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))

        title='_'+attr+'_frame_list_union_pivot.parquet'
        contacts_pivot_full_list.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        #title='_'+attr+'_frame_list_union_pivot.pkl'
        #contacts_pivot_full_list.to_pickle(f'{data_path}/{label}{title}')
        

        #list[len(C_k(i,j))] ---> distribution of percistance
        contacts_pivot_list_perc = dt_cont_map_all_sym.pivot_table(
        index=attr+"_a",
        columns=attr+"_b",
        values='percistance_times',
        aggfunc=make_list,
        observed=False
        )
        contacts_pivot_full_list_perc = contacts_pivot_list_perc.reindex(index=attr_val, columns=attr_val,fill_value=0) 
        contacts_pivot_full_list_perc=contacts_pivot_full_list_perc.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))
        

        title='_'+attr+'_list_percistance_distribution_pivot.parquet'
        contacts_pivot_full_list_perc.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        
        #title='_'+attr+'_list_percistance_distribution_pivot.pkl'
        #contacts_pivot_full_list_perc.to_pickle(f'{data_path}/{label}{title}')

        if bool_bonds:
            #special bond_percistance time 
            for key in ['cation_pi_t','pi_stacking_t','salt_bridge_t']:
    
                key_c=key+'_pt'
    
                #dt_cont_map_all_sym[key_c]=dt_cont_map_all_sym[key].apply(lambda x: [len(sub) for sub in split_consecutive_intervals(x)] if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))
    
                dt_cont_map_all_sym[key_c] = dt_cont_map_all_sym[key].apply(lambda x: [len(sub) for sub in split_consecutive_intervals(x)])
            
                #pivot
                #attr_a='res_type_a'
                #attr_b='res_type_b'
                contacts_pivot_list_perc_bond = dt_cont_map_all_sym.pivot_table(
                        index=attr+"_a",
                        columns=attr+"_b",
                        values=key_c,
                        aggfunc=make_list,
                        observed=False
                        )
                
                contacts_pivot_full_list_perc_bond = contacts_pivot_list_perc_bond.reindex(
                    index=attr_val, columns=attr_val, fill_value=0).fillna(0)
            
                contacts_pivot_full_list_perc_bond=contacts_pivot_full_list_perc_bond.applymap(lambda x: x if isinstance(x, list) else ([] if x is None or (isinstance(x, float) and pd.isna(x)) else []))
    
                title='_'+attr+'_'+key+'_list_percistance_distribution_pivot.parquet'
                contacts_pivot_full_list_perc_bond.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
    
    
                #frequency
            
                #sums
                col=key
                # Create a new column with the suffix '_sum' (e.g., 'cation_pi_t_sum')
                sum_col = f"{col}_sum"
                
                # Calculate the length of the list in each cell; 0 if not a list/array/string
                dt_cont_map_all_sym[sum_col] = dt_cont_map_all_sym[col].apply(
                    lambda x: len(x) if isinstance(x, (list, tuple, str, np.ndarray)) else 0
                )
    
                #manual sol:
                pivot_name = f"contacts_pivot_{col}_freq"
                
                # Create the base pivot table
                base_pivot = dt_cont_map_all_sym.pivot_table(
                    index=attr+"_a",
                    columns=attr+"_b",
                    values=sum_col,
                    aggfunc='sum', # Use 'sum' from numpy/pandas
                    observed=False
                )
                
                base_pivot_full_time = base_pivot.reindex(index=attr_val, columns=attr_val,fill_value=0)
    
                
                title='_'+attr+'_'+key+'_time_count_pivot.parquet'
                base_pivot_full_time.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")


        #avrg(list[len(C_k(i,j))]) ---> average from distribution of percistance
        contacts_pivot_full_list_perc_avrg=contacts_pivot_full_list_perc.applymap(
    lambda x: np.mean(x) if isinstance(x, list) else x).fillna(0)

        title='_'+attr+'_avrg_percistance_distribution_pivot.parquet'
        contacts_pivot_full_list_perc_avrg.to_parquet(f'{data_path}/{label}{title}',engine="pyarrow")
        return 0


def main(args):


        
    path = args.path
    path_tools = args.path_tools
    output_path = args.output
    contact_path=output_path
    label=args.label
    results_name= args.results_name
    bool_pwi= bool(strtobool(args.bool_pwi))
    bool_bonds=bool(strtobool(args.bool_bonds))
    
    #################### load data ##############################
    #############################################################
    
    data_path=contact_path+"/"+results_name
    title='_contacts_only_t.parquet' 
    
    if bool_bonds:
         dt_cont_map_all_sym=pd.read_parquet(f'{data_path}/{label}{title}', columns=['res_a', 'res_b', 'frame','frame_lb','frame_d2t','len_contacts','mol_ind_a', 'res_dom_a', 'res_type_a', 'struc_id_a', 'prot_a','res_org_a', 'mol_ind_b', 'res_dom_b', 'res_type_b', 'struc_id_b','prot_b', 'res_org_b','percistance_times','avg_persistence', 'median_persistence', 'cation_pi_t', 'pi_stacking_t', 'hbond_t', 'salt_bridge_t'])
    else:
        dt_cont_map_all_sym=pd.read_parquet(f'{data_path}/{label}{title}', columns=['res_a', 'res_b', 'frame','frame_lb','frame_d2t','len_contacts','mol_ind_a', 'res_dom_a', 'res_type_a', 'struc_id_a', 'prot_a','res_org_a', 'mol_ind_b', 'res_dom_b', 'res_type_b', 'struc_id_b','prot_b', 'res_org_b','percistance_times','avg_persistence', 'median_persistence'])


    #################### load metadata ##########################
    #############################################################
    
    title = label+'_contacts_res_meta'
    res_df=pd.read_parquet(contact_path + '/' + title + '.parquet')
    

    #################### evaluate data ##########################
    #############################################################
    
    #attributes to evaluate
    l_attributes=['res_type','res_org','prot']


    #################### single attributes ######################
    #############################################################

    l_attr_values=[]
    for attr in l_attributes:

        #calc unnique values of the attribute (e.g. all residue type pairs, or all combination of protein types
        if attr=='res_type':
            l_attr_values.append(amino_acids)
        else:
            l_attr_values.append(sorted(res_df[attr].dropna().unique()))
        
    #run
    for attr, attr_val in zip(l_attributes,l_attr_values):
        
        analyse_attribute(attr, attr_val,dt_cont_map_all_sym,data_path,args)


if __name__ == "__main__":
    
    args = get_args()
    main(args)