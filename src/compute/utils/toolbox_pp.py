import re
import os
import shutil
import glob
import itertools
import numpy as np
import pickle
import pandas as pd
import MDAnalysis as mda
from Bio.Seq import Seq, MutableSeq
import mdtraj as md

#find concecutive index sets in list of overall indicies
def intervals_extract(iterable):
    #iterable = sorted(iterable)
    for key, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [item[1] for item in group]

#find all occurances of seq_cut in seq_multi
def find_all_occurrences(seq_multi, seq_cut):
    start = 0
    occurrences = []
    
    while True:
        start = seq_multi.find(seq_cut, start)
        if start == -1:
            break
        occurrences.append(start)
        start += 1  # Move past the last found position to find subsequent occurrences
    return occurrences

#multimer into units - by sequence comparison
def calc_multimer_components(df_domains, prot):
    #path_input_a="/home/lubaltz/code/SFB1551/r08/GMX_models/inputs/PEI1_PEI2"
    #monomers
    namA=df_domains["prot"][0]#"PEI1"
    namB=df_domains["prot"][1]#"PEI2"
    #full proteins & cuts
    A_top=df_domains[df_domains["prot"]==namA]#path_input_a+"/"+"PEI1_07d79.pdb" #PEI1C
    B_top=df_domains[df_domains["prot"]==namB]#path_input_a+"/"+"PEI2_cd323.pdb" #PEI2C
    u_full_A=mda.Universe(A_top['file'][0])
    u_full_B=mda.Universe(B_top['file'][1])
    cutA=u_full_A.residues.resnames[A_top['min'][0]-1:A_top['max'][0]]
    cutB=u_full_B.residues.resnames[B_top['min'][1]-1:B_top['max'][1]]

    #extract sequences
    l_cutA=[three_to_one(item) for item in cutA]
    str_cutA=''.join(l_cutA)
    l_cutB=[three_to_one(item) for item in cutB]
    str_cutB=''.join(l_cutB)
    
    multi=prot.residues.resnames
    l_multi=[three_to_one(item) for item in multi]
    str_multi=''.join(l_multi)

    print("str_multi", str_multi)

    print("str_A", str_cutA)
    print("str_A", str_cutB)

    #find sequences
    seq_cutA = Seq(str_cutA)
    seq_cutB = Seq(str_cutB)
    seq_multi= Seq(str_multi)

    ind_a=find_all_occurrences(seq_multi, seq_cutA)
    ind_b=find_all_occurrences(seq_multi, seq_cutB)

    
    l_ind = sorted(ind_a + ind_b)
    l_nam = [
        namA if item in ind_a else namB
        for item in l_ind
    ]
    print("l_ind", l_ind)
    return l_nam


#expand domain df with multimer data
def expand_df_domains(molA, molB, prot_org,top_mix,names, df_domains, label_new,bookkeeping_path, multimers=True):
    multi_A=mda.Universe(str(molA)+'.pdb')
    multi_B=mda.Universe(str(molB)+'.pdb')

    
    
    l_comp_A=calc_multimer_components(df_domains, multi_A)
    l_comp_B=calc_multimer_components(df_domains, multi_B)

    print(l_comp_A, "comp_A")
    print(l_comp_B, "comp_B")
    dict_A={'prot':"-".join(l_comp_A),'cut_na':len(multi_A.atoms) ,'file':str(molA)+'.pdb', 'multimer':l_comp_A}
    dict_B={'prot':"-".join(l_comp_B),'cut_na':len(multi_B.atoms) ,'file':str(molB)+'.pdb', 'multimer':l_comp_B}
    print(dict_A, "dict_A")
    print(dict_A, "dict_B")
    if dict_A==dict_B:
        df_dict = pd.DataFrame([dict_A])
    else:
        df_dict = pd.DataFrame([dict_A, dict_B])

    if multimers:
        df_domains_new = pd.concat([df_domains, df_dict], ignore_index=True)
        df_domains_new[['orig_na','cut_na','min','max' ,'BTB_min' ,'BTB_max' ,'BACK_min','BACK_max']]=df_domains_new[['orig_na','cut_na','min','max' ,'BTB_min' ,'BTB_max' ,'BACK_min','BACK_max']].fillna(0)
        df_domains_new[['orig_na','cut_na','min','max' ,'BTB_min' ,'BTB_max' ,'BACK_min','BACK_max']] = df_domains_new[['orig_na','cut_na','min','max' ,'BTB_min' ,'BTB_max' ,'BACK_min','BACK_max']].astype(int)
    
        df_domains_new=df_domains_new.drop_duplicates(subset='prot', keep='first').reset_index(drop=True)
    #print("domains_new",df_domains_new)
        sys_dict_new=setup_sys_dict_ids(prot_org, top_mix, names, df_domains=df_domains_new) #for multimers
    
    else:
        df_domains_new=df_domains
        sys_dict_new=setup_sys_dict_ids(prot_org, top_mix, names)
        
    #
    
    #sys_dict_new=setup_sys_dict_ids(prot_org, top_mix, names)
    #print("sys_dict_new",sys_dict_new)
    l_dict, sys_dict_end=calc_sys_domains(sys_dict_new, df_domains_new)
    #print("sys_dict_end",sys_dict_end)
    #print("l_dict",l_dict)
    #domains_sys_new=pd.DataFrame(l_dict)

    return l_dict, sys_dict_end


#from topology file extract a list of system compounts 
def extract_compound_list(filename):
    compound_dict = {}
    with open(filename, 'r') as file:
        lines = file.readlines()
        compound_section = False
        for line in lines:
            if line.strip().startswith('[ molecules ]'):
                compound_section = True
                continue
            if compound_section and re.match(r'\s*\w+\s+\d+', line):
                compound, count = line.strip().split()
                compound_dict[compound] = int(count)
    return compound_dict

def gather_xtcs(path, key='*', bool_pp=True, name="md_run"):
    if os.path.isfile(os.path.join(path, "md_run.part0001.xtc")): 
        pass
    elif bool_pp:
        shutil.copy(os.path.join(path, name+".xtc"), os.path.join(path, name+".part0001.xtc"))
    xtc_files = glob.glob(os.path.join(path, key+"part*.xtc"))
    xtc_sorted=sorted(xtc_files)
    return xtc_sorted


def setup_sys_dict_ids(prot_org, top_sys, names, df_domains=0):

    #protein ids
    prot_mols=list(intervals_extract(prot_org.residues.resids)) #mols in total
    print("prot mols", prot_mols)
    print("there exits", len(prot_mols), "protein chains")
    n_mols=[len(item) for item in prot_mols]
    
    #protein names
    compound_dict = extract_compound_list(top_sys)
    #print("compounds in .top",compound_dict)
    #protein_names_org = [key.split('_')[-1] for key, value in compound_dict.items() if key.startswith('Protein') for _ in range(value)]
    protein_names_org = [key.split('_')[-1] for key, value in compound_dict.items() if key.startswith(('Protein', 'Molecule')) for _ in range(value)]
    protein_names = [names.get(item, item) for item in protein_names_org]
    print("names",protein_names)

    if type(df_domains)!=int: #.empty!=True:
        #print("shift")
        
        n_mols=[]
        sys_protein_names=[]
        monomer_names=[]
        orig_min=[]
        prot_mols_all=[]
        m_mols=0
        l_mols=[]
        #all molecules in the system (can be multimers or monomers)
        #print(protein_names, prot_mols)
        for mol_name, mol_ids in zip(protein_names, prot_mols):
            #sys_protein_names.append(mol_name) 
            print("mol_name", mol_name)
            print("domains", df_domains)
            multi=df_domains[df_domains['prot']==mol_name]['multimer'].values[0]
            #print(multi)
            if isinstance(multi, list): #expand to monomers
                if len(multi)>1 :
                    for mon in multi:
                        #print("hi")
                        n_mols.append(len(mol_ids))
                        sys_protein_names.append(mol_name)
                        monomer_names.append(mon) 
                        orig_min=int(df_domains[df_domains['prot']==mon]['min'])
                        prot_mols_all.append(mol_ids+(orig_min-1)*np.ones(len(mol_ids), dtype='int'))
                        l_mols.append(m_mols)
            else:
                n_mols.append(len(mol_ids))
                sys_protein_names.append(mol_name)
                #print(df_domains[df_domains['prot']==mol_name])
                orig_min=int(df_domains[df_domains['prot']==mol_name]['min']) #minimum original resid
                prot_mols_all.append(mol_ids+(orig_min-1)*np.ones(len(mol_ids), dtype='int')) #for min, prot in zip(orig_min,prot_mols)]
                monomer_names.append(mol_name)
                l_mols.append(m_mols)
                
            #shifted ids according to cut given by shif_min if not given by 
            #orig_min=[int(df_domains[df_domains['prot']==item]['min']) for item in protein_names]
            #prot_mols=[prot+(min-1)*np.ones(len(prot), dtype='int') for min, prot in zip(orig_min,prot_mols)]
            m_mols=m_mols+1
    #else:
    #    n_mols=[len(item) for item in prot_mols]
    #    #prot_mols_all=
    #    #monomer_names=
    #    #sys_protein_names=
    #    #l_mols=
    dict={'n_res':n_mols, 'l_org_mol_ids':prot_mols_all, 'type_mols':monomer_names, "sys_mols":sys_protein_names, 'ind_mols':l_mols}
    return dict

def calc_sys_index(n,type, curr_end, domains, offset,mol_id):
    d_line={}
    mol_diff=offset
    mol_min=domains[domains['prot']==type]['min'].iloc[0]
    mol_max=domains[domains['prot']==type]['max'].iloc[0]
    BTB_min=domains[domains['prot']==type]['BTB_min'].iloc[0]
    BTB_max=domains[domains['prot']==type]['BTB_max'].iloc[0]
    BACK_min=domains[domains['prot']==type]['BACK_min'].iloc[0]
    BACK_max=domains[domains['prot']==type]['BACK_max'].iloc[0]

    mol_min_t=mol_min -mol_diff + curr_end
    mol_max_t=mol_max -mol_diff + curr_end
    BTB_min_t=BTB_min -mol_diff + curr_end 
    BTB_max_t=BTB_max -mol_diff + curr_end
    BACK_min_t=BACK_min -mol_diff + curr_end
    BACK_max_t=BACK_max -mol_diff+ curr_end

    d_line={'prot':type,'min': mol_min_t,'max':mol_max_t,'BTB_min':BTB_min_t,'BTB_max':BTB_max_t,'BACK_min':BACK_min_t,'BACK_max': BACK_max_t, 'ind_mol':mol_id}
    
    str_mol=f'resid {mol_min -mol_diff +curr_end}-{mol_max -mol_diff +curr_end}'
    str_btb_mol= f'resid {BTB_min -mol_diff + curr_end }-{BTB_max -mol_diff + curr_end}'
    str_back_mol= f'resid {BACK_min -mol_diff + curr_end }-{BACK_max-mol_diff+ curr_end}'
    curr_end=curr_end+mol_max-mol_min+1

    return curr_end,str_mol, str_btb_mol, str_back_mol, d_line


def calc_sys_domains(sys_dict, df_domains):
    #domain information
    prot_mols=sys_dict['l_org_mol_ids'] #protein index from input files 1:max_a, 1:max_b
    #print("there exits", sys_dict['n_res'], "protein chains residues")
    a=prot_mols[0] #orginal res ids from full sequence
    b=prot_mols[1] #orginal res ids from full sequence
    diff_a=a[0]-1 #starting position a in full A sequence
    diff_b=b[0]-1 #starting index b in full B sequence

    end=0
    l_mol=[]
    l_btb=[]
    l_back=[]
    l_dict=[]
    l_mol_id=[]
    #u_mol=[]
    #u_btb=[]
    #u_back=[]
    for i, item in enumerate(sys_dict['type_mols']):
        print (i,item)
        mol_id=sys_dict['ind_mols'][i]
        delta=sys_dict['l_org_mol_ids'][i][0]-1
        end, mol, btb, back, d_line=calc_sys_index(i,item, end, df_domains, delta, mol_id)
        #print(end, mol, btb, back)
        l_dict.append(d_line)
        l_mol.append(mol)
        l_btb.append(btb)
        l_back.append(back)
        #l_mol_id.append()
        #u_mol.append(u.select_atoms(mol))
        #u_btb.append(u.select_atoms(btb))
        #u_back.append(u.select_atoms(back))
    
    #sys_dict['u_mol']=u_mol
    sys_dict['mol']=l_mol
    sys_dict['btb']=l_btb
    sys_dict['back']=l_back
    #sys_dict['ind_mol']=l_mol_id

    return l_dict, sys_dict

def ps_to(ps, to='ns'):
    
    if to=='ns':
        return ps/1e+3
    elif to=='mic':
        return ps/1e+6 
    elif to=='mil':
        return ps/1e+9
    elif to=='sec':
        return ps/1e+12
    else:
        print('unknown unit')
        return 0

def create_bash_file(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel              	        # Partition/Queue name
#SBATCH -C skylake #broadwell #skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N {1:01}                     	        #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc             	   # Specify allocation to charge against
##SBATCH --mem=100000
##SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))

    # Append each command line by line
    if isinstance(command, list):
        for cmd in command:
            content += cmd.strip() + "\n"
    else:
        content += command.strip() + "\n" # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content




def create_bash_file_memory(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel              	        # Partition/Queue name
#SBATCH -C broadwell #broadwell #skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N {1:01}                     	        #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc             	   # Specify allocation to charge against
#SBATCH --mem=246000
#SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))

    # Append each command line by line
    if isinstance(command, list):
        for cmd in command:
            content += cmd.strip() + "\n"
    else:
        content += command.strip() + "\n" # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content


def create_bash_file_high_memory2(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p bigmem           # Partition/Queue name 'bigmem'
#SBATCH -C broadwell          # select either 'broadwell' or 'skylake'
#SBATCH -N {1:01}           #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00     #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=800000
#SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))

    # Append each command line by line
    if isinstance(command, list):
        for cmd in command:
            content += cmd.strip() + "\n"
    else:
        content += command.strip() + "\n" # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content





def create_bash_file_high_memory(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p bigmem           # Partition/Queue name 'bigmem'
#SBATCH -C broadwell          # select either 'broadwell' or 'skylake'
#SBATCH -N {1:01}           #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00     #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=0
#SBATCH --exclusive

source /home/lubaltz/p_3_10/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))

    # Append each command line by line
    if isinstance(command, list):
        for cmd in command:
            content += cmd.strip() + "\n"
    else:
        content += command.strip() + "\n" # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content

def create_bash_file_medium_memory(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p bigmem           # Partition/Queue name 'bigmem'
#SBATCH -C skylake          # select either 'broadwell' or 'skylake'
#SBATCH -N {1:01}           #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00     #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A m2_komet331hpc   # Specify allocation to charge against
#SBATCH --mem=0
#SBATCH --exclusive


source /home/lubaltz/p_3_10/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))

    # Append each command line by line
    if isinstance(command, list):
        for cmd in command:
            content += cmd.strip() + "\n"
    else:
        content += command.strip() + "\n" # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content



def create_bash_file_NHR(command, hours=24, nodes=4, name='bla'):
    content = """#!/bin/bash
#-----------------------------------------------------------------
# {0}
# 
#-----------------------------------------------------------------

#SBATCH -J "{0}"         	#todo Job name #todo
#SBATCH -o "{0}".%j.out  	#todo Specify stdout output file (%j expands to jobId)
#SBATCH -p parallel              	        # Partition/Queue name
#SBATCH -N {1:01}                     	        #todo 8 Total number of nodes requested (64 cores/node)
#SBATCH -t {2:02}:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
#SBATCH -A nhr-dyndisphase             	   # Specify allocation to charge against

source /home/lubaltz/nhr_tools/bin/activate

# Commands
""".format(name,nodes,hours) # print('{:01}:{:01}:{:01}'.format(int(3), int(4), int(5)))
    
    content += command.strip()  # Ensure the command doesn't start with a newline
    content += """

echo byebye bash
""" 
    return content


def format_time(days, hours, minutes, seconds):
    # Ensure all time components are two digits where necessary
    return f"{days}-{hours:02}:{minutes:02}:{seconds:02}"

#def create_bash_file(output_path, comand):
#    #VMD Labeling
#    #write domains to VMD file tcl
#    mols, btb, back, betas=calc_vmd_labels(domains_sys)
#    
#    
#    content = """
##    #!/bin/bash

    #-----------------------------------------------------------------
    # submitting papermill
    # 
    #-----------------------------------------------------------------
    
    #SBATCH -J "submitted_papermill"         	#todo Job name #todo
    #SBATCH -o "submitted_papermill".%j.out  	#todo Specify stdout output file (%j expands to jobId)
    #SBATCH -p parallel              	        # Partition/Queue name
    #SBATCH -C skylake          		        # select either 'broadwell' or 'skylake'
    #SBATCH -N 4                     	        #todo 8 Total number of nodes requested (64 cores/node)
    #SBATCH -t 24:00:00              	        #todo  Run time (hh:mm:ss) - 0.5 hours
    #SBATCH -A m2_komet331hpc             	   # Specify allocation to charge against

#    source /home/lubaltz/p_3_10/bin/activate
    
    # comands
#    """.format(top)
#    content +=comand
#    content += """
    
#    byebye bash
#   """ 
#    return content

def three_to_one(resname):
    # Dictionary mapping three-letter codes to one-letter codes
    three_to_one_dict = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    one_letter_code = three_to_one_dict.get(resname, 'X')  # 'X' for unknown or non-standard residues

    return one_letter_code

def find_and_sort_min_max_columns(df):
    """
    Find all column name groups that have both *_min and *_max variants,
    and sort them by the ascending order of the minimum value
    found in each *_min column.
    
    Parameters
    ----------
    df : pandas.DataFrame
    
    Returns
    -------
    list
        List of base names (without _min/_max), sorted by the column minima.
    """
    # Find base names that have both *_min and *_max columns
    cols = df.columns
    min_cols = {c[:-4] for c in cols if c.endswith('_min')}
    max_cols = {c[:-4] for c in cols if c.endswith('_max')}
    common = sorted(list(set(min_cols & max_cols)))

    # Get the true minimum of each *_min column
    sort_values = {base: df[f"{base}_min"].min() for base in common}

    # Sort by ascending order of the *_min column’s minimum value
    sorted_bases = [k for k, _ in sorted(sort_values.items(), key=lambda x: x[1])]
    
    return sorted_bases

def calc_r0(struc, radius,df_sys_domains_file,NATIVE_CUTOFF=6): #from dimer

    df_sys_domains=pd.read_parquet(df_sys_domains_file)


    prots=list(set(df_sys_domains["prot"]))
    #reference dimer
    u_ref = mda.Universe(struc)
    
    #identify
    dat_A=df_sys_domains.iloc[0]
    dat_B=df_sys_domains.iloc[1]
    multi_A=u_ref.select_atoms(f"resid {dat_A['min']}-{dat_A['max']}").atoms.ids-1
    multi_B=u_ref.select_atoms(f"resid {dat_B['min']}-{dat_B['max']}").atoms.ids-1
    min_multi_A=multi_A.min()
    min_multi_B=multi_B.min()

    
    #load inter chain pairs
    combinations_nat = np.array(list(itertools.product(multi_A, multi_B)), dtype=int)

    #check distances for inter chain pairs
    ref_frame=md.load_pdb(struc)
    ref_pairs_distances = md.compute_distances(ref_frame, combinations_nat,periodic=True)[0]
    native_contacts = combinations_nat[ref_pairs_distances < NATIVE_CUTOFF]

    #load inter chain pairs
    r0 = md.compute_distances(ref_frame, native_contacts,periodic=True)

    #save atom indicies (wrt. individual chains)
    nat_A=native_contacts[:,0]
    nat_B=native_contacts[:,1]
    
    #native index w.r.t reference monomer structure
    ind_nat_A=nat_A-min_multi_A
    ind_nat_B=nat_B-min_multi_B

    return ind_nat_A, ind_nat_B, r0, prots
