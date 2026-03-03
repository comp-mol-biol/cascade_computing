import re
import os
import shutil
import glob
import itertools
import numpy as np
import MDAnalysis as mda


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

def one_to_three(one_letter):
    # Dictionary mapping one-letter codes to three-letter codes
    one_to_three_dict = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
        'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
        'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
        'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
        'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    
    return one_to_three_dict.get(one_letter, 'UNK')

def guess_iter_steps (ns_per_day, iter_days, dt):

    #per day
    ps_per_day=ns_per_day*1000
    steps_per_day=ps_per_day/dt

    return iter_days*steps_per_day

def calc_rho_p(mols,n_mols, lx, ly, lz, d):
    all_atm=0
    v=(lx+d)*(ly+d)*(lz+d)
    for mol, n in zip(mols, n_mols):
        u_mol=mda.Universe(mol)
        n_atm=n*len(u_mol.atoms)
        all_atm+=n_atm
        
    return all_atm/v



def check_sim(md_file): 
    step_value=0
    with open(md_file, 'r') as file:
        #lines = file.readlines()
        lines = [line.strip() for line in file if line.strip()]
    if lines and lines[-1].startswith('Finished mdrun'):
        print(lines[-1])
        k=0
        while k < len(lines):
            if lines[-1-k].startswith('Writing checkpoint, step '):
                match = re.search(r'step\s+(\d+)', lines[-1-k])
                if match:
                    step_value = match.group(1)
                    #print(f"Extracted step value: {step_value}")
                break
            k=k+1
            
        return True, step_value
    elif lines and not lines[-1].startswith('Finished mdrun'):
        print(lines[-1])
        k=0
        while k < len(lines):
            if lines[-1-k].startswith('Writing checkpoint, step '):
                match = re.search(r'step\s+(\d+)', lines[-1-k])
                if match:
                    step_value = match.group(1)
                    #print(f"Extracted step value_no done: {step_value}")
                break
            k=k+1
            #print(step_value)
        return False, step_value
    else:
        return False, 0


def check_run(md_file, l=-1, key='Finished mdrun'): 
    with open(md_file, 'r') as file:
        lines = file.readlines()
    if lines and lines[l].startswith(key):
        return True
    else:
        return False

#simplify copying
def copy_and_rename_file(src_directory, dest_directory, label, pattern, all=False):
    # Search for files starting with 'Q_x' in the source directory
    search_pattern = os.path.join(src_directory, pattern)
    files = glob.glob(search_pattern)

    if not files:
        print("No file found.")
        return None

    # Get the first matching file
    src_file = files[0]

    if all:
        for src_file in files:
            # Extract the file name without the directory path
            file_name = os.path.basename(src_file)
        
            # Extract the file name without extension and the extension
            file_name_without_ext, file_ext = os.path.splitext(file_name)
        
            # Create the new file name with the label
            new_file_name = f"{label}{file_name_without_ext}{file_ext}"
        
            # Construct the full path of the destination file
            dest_file = os.path.join(dest_directory, new_file_name)
        
            # Copy the file to the new directory with the new name
            shutil.copy(src_file, dest_file)
            print(f"File {file_name} copied and renamed to {new_file_name} in {dest_directory}")       

    else:

        # Extract the file name without the directory path
        file_name = os.path.basename(src_file)
    
        # Extract the file name without extension and the extension
        file_name_without_ext, file_ext = os.path.splitext(file_name)
    
        # Create the new file name with the label
        new_file_name = f"{label}{file_name_without_ext}{file_ext}"
    
        # Construct the full path of the destination file
        dest_file = os.path.join(dest_directory, new_file_name)
    
        # Copy the file to the new directory with the new name
        shutil.copy(src_file, dest_file)
        print(f"File {file_name} copied and renamed to {new_file_name} in {dest_directory}")

    return dest_file