import os
import re
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

amino_acids = [
    "ALA","ARG","ASN","ASP","CYS",
    "GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO",
    "SER","THR","TRP","TYR","VAL"
]

def query_slice(query, subquery,matT):
    q1=query[0]
    q2=query[1]
    s1=subquery[:,0]
    s2=subquery[:,1]
    mat=matT.loc[s1, s2]
    return mat

def calc_frequency(time_contacts_prot_res_type,key_a='', key_b=''):
    contact_frequency=time_contacts_prot_res_type.copy()

    if isinstance(time_contacts_prot_res_type.index, pd.MultiIndex):
    
        contact_frequency.loc[(key_a, key_b)]=time_contacts_prot_res_type.loc[(key_a, key_b)].to_numpy()/float(time_contacts_prot_res_type.to_numpy().sum())

    else:
        contact_frequency=time_contacts_prot_res_type/float(time_contacts_prot_res_type.to_numpy().sum())

    return contact_frequency


def normalize_frequency(time_contacts_prot_res_type,int_res_types,key_a='', key_b=''):

    contact_frequency_normalized=time_contacts_prot_res_type.copy()

    if isinstance(time_contacts_prot_res_type.index, pd.MultiIndex):
        #normalized count
        numerator=contact_frequency_normalized.loc[(key_a, key_b)]
        denominator=int_res_types
        val = numerator.div(denominator.reindex(numerator.index).replace(0, np.nan)).fillna(0)/time_contacts_prot_res_type.to_numpy().sum()
        contact_frequency_normalized.loc[(key_a, key_b)]=val.to_numpy()

    else: 
        numerator=contact_frequency_normalized
        denominator=int_res_types
        val = numerator.div(denominator.reindex(numerator.index).replace(0, np.nan)).fillna(0)/time_contacts_prot_res_type.to_numpy().sum()
        contact_frequency_normalized=val#.to_numpy()

    return contact_frequency_normalized


def plot_df(df,xlabel=" ",ylabel=" ",label="", file='./plot'):
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    heatmap=sbn.heatmap(df,xticklabels=True,yticklabels=True,square=False,cbar=True,cbar_kws={"shrink": 0.75},cmap="Spectral_r",annot=False)
    heatmap.invert_yaxis()
    colorbar = heatmap.collections[0].colorbar  # Get the colorbar from the heatmap
    colorbar.set_label(label)
    colorbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.4f'))
    fig.tight_layout()
    ax.set_xlabel(xlabel)#, fontsize=32)
    ax.set_ylabel(ylabel)#, fontsize=32)
    # Increase tick label font size
    ax.tick_params(axis='x',  rotation=90)
    ax.tick_params(axis='y')
    print(file)
    fig.savefig(file, dpi=100, bbox_inches="tight")


#function: - evaluate aa table
def eval_aa_table(res_df, name='MUT16_chain'):

    grouped_p = (
        res_df.groupby(['struc_id','res_type'])
        .size()
        .reset_index(name='count')
    )
    pivot_df_chains = grouped_p.pivot_table(
        index='res_type',
        columns='struc_id',
        values='count',
        fill_value=0,  # Use 0 for struc/res_type combos that don't exist
        observed=False
    )
    
    grouped_p = (
        res_df.groupby(['prot','res_type'])
        .size()
        .reset_index(name='count')
    )
    res_types_df = grouped_p.pivot_table(
        index='res_type',
        columns='prot',
        values='count',
        fill_value=0,  #combos that don't exist,
        observed=False
    )
    
    
    res_types_df['total'] = res_types_df.sum(axis=1)
    res_types_df[name] =pivot_df_chains.iloc[:, 0]
    res_types_df = res_types_df.reindex(amino_acids, fill_value=0)

    return res_types_df