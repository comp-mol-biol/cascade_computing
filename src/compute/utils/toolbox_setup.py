import sys
import MDAnalysis as mda
import pandas as pd
import re
import glob
import os
from .toolbox_sim import three_to_one


def gather_xtcs(path, key='*', bool_pp=True, name="md_run", part=True):
    if os.path.isfile(os.path.join(path, "md_run.part0001.xtc")): 
        pass
    elif bool_pp:
        shutil.copy(os.path.join(path, name+".xtc"), os.path.join(path, name+".part0001.xtc"))
    if part:
        xtc_files = glob.glob(os.path.join(path, key+"part*.xtc"))
    else:
        xtc_files = glob.glob(os.path.join(path, key+"*.xtc"))
    xtc_sorted=sorted(xtc_files)
    return xtc_sorted

def add_sequences(df_domains):

    #full sequence
    df_domains["full_seq"] = df_domains["file"].apply(
        lambda f: "".join([three_to_one(res.resname) for res in mda.Universe(f).residues]))
    
    #cut sequence
    df_domains["cut_seq"] = df_domains.apply(
        lambda row: "".join([
            three_to_one(res.resname) 
            for res in mda.Universe(row["file"]).residues[row["min"]-1 : row["max"]]
        ]), 
        axis=1  # Apply to rows, not columns
    )
    return df_domains


def find_sequence_pos(prot_sequence, sequence):
    # Find all matches of the sequence in the prot_sequence
    matches = re.finditer(f'(?={re.escape(sequence)})', prot_sequence)
    # Collect start and end positions for each match
    positions = [(match.start(), match.start() + len(sequence)-1) for match in matches]
    return positions


#compare with monomers (what chain type is a segment?)
def make_mon_df(top,t,df_domains):
    u=mda.Universe(top)
    sys=pd.DataFrame
    df_mons = pd.DataFrame(columns=["prot", "struc_id"])
    
    for l,item in enumerate(u.segments):
        print("seg",item)
        #seq = [three_to_one(res.resname) for res in item.residues if not set(res.atoms.names).issubset({"CA"})]#if res.atoms.names[0] != "CA"]
        seq=[three_to_one(res.resname) for res in item.residues if not set(res.atoms.names).issubset({"CA"})]
        str_seq= ''.join(seq)
        print("l seg",len(seq))
    
        #find 
        for k, row in df_domains.iterrows():
            print("monomer_unit", k, row["cut_seq"])
            print("chain", k, str_seq)
            l_mon=find_sequence_pos(str_seq,row["cut_seq"])
            for j in l_mon:
                df_mons.loc[len(df_mons)] = [row["prot"], t]
    print(df_mons)
    return df_mons


#get subunits from input file (pdb - inital unit simulations)

def make_sys_mon_df(filenames,df_domains):
    sys_mon_df=pd.DataFrame(columns=["prot", "struc_id"])
    for t, file in enumerate(filenames):
        print("inside",t,file)
        #top=path_input+"/"+file+".pdb"
        tmp_df_mons=make_mon_df(file[0],t,df_domains)
        tmp_df_mons["unit"]=file[1]
        sys_mon_df=pd.concat([sys_mon_df, tmp_df_mons])

    return sys_mon_df


# create sys dicts - wrt monomers provided in df_domains["prot"]

def make_sys_mon_df_max(sys_mon_df,df_domains):
    max_prev = 0  # Initialize max_prev before the loop
    sys_mon_df.reset_index(drop=True, inplace=True)
    for i, row in sys_mon_df.iterrows():
        # Find the match for the current row
        match = df_domains[df_domains["prot"] == row["prot"]]
        
        # Calculate the min and max index values
        min_ind = match["min"].values[0] - 1     #1 ->0
        max_ind = match["max"].values[0] - 1     #5 ->4maximum not included
        len_mon = max_ind - min_ind + 1          #5
        
        # Update the "min" and "max" residue ids (not indicies) values in the DataFrame using .loc[]
        sys_mon_df.loc[i, "min"] = max_prev + 1  
        sys_mon_df.loc[i, "max"] = max_prev + len_mon
    
        #sys_mon_df.loc[i, "BTB_min"] = max_prev
        
        # Update max_prev to be the "max" value of the current row
        max_prev = sys_mon_df.loc[i, "max"]
    return sys_mon_df


# create sys dicts - wrt to domains provided in 


def make_sys_mon_df_max_domains(sys_mon_df,df_domains,domains):
    #iterate over monomers
    sys_mon_df.reset_index(drop=True, inplace=True)
    for i, row in sys_mon_df.iterrows():
    
        #find domains in df_domais for every monomer type
        match = df_domains[df_domains["prot"] == row["prot"]]
    
        for dom in domains:
            dom_min= dom+"_min"
            dom_max= dom+"_max"
    
            try:
                
                sys_mon_df.loc[i, dom_min]=row["min"]+match[dom_min].values[0] - match["min"].values[0]
                sys_mon_df.loc[i, dom_max]=row["min"]+match[dom_max].values[0] - match["min"].values[0]
            except:
                print("no domain ", dom, "in ",row["prot"]) 
        selected_cols = [col for col in sys_mon_df.columns if col not in ["prot", "unit"]]
        for col in sys_mon_df[selected_cols]:
            sys_mon_df[col]= sys_mon_df[col].fillna(0).astype(int)
    return sys_mon_df

def calc_sys_domain(molA, molB, df_domains, n_mols):
    df_domains=add_sequences(df_domains)
    geo_mols=[molA, molB]
    labels=["A", "B"]
    files=[[item, lab] for count, item, lab in zip(n_mols, geo_mols, labels) for _ in range(count)]

    df_sys_domains=make_sys_mon_df(files,df_domains)
    df_sys_domains=make_sys_mon_df_max(df_sys_domains,df_domains)
    df_sys_domains=make_sys_mon_df_max_domains(df_sys_domains,df_domains,domains)
    
    return df_sys_domains