"""
MDAnalysis Trajectory Processing Script.

Description
-----------
This script performs post-processing on Gromacs molecular dynamics trajectories.
It operates sequentially to load topology and trajectory data, filter for
protein and ion atoms, assign chain identifiers, and apply coordinate
transformations (centering and wrapping).

The script executes the following steps:
1.  **Configuration**: Sets paths and parameters.
2.  **Metadata Loading**: Reads topology (.tpr) and trajectory lists (.pkl).
3.  **Initial Extraction**: Saves the initial protein/ion structure as a PDB.
4.  **Trajectory Loading**: Loads the full trajectory (out-of-core).
5.  **Atom Selection**: Filters for protein and specific ions, excluding virtual sites.
6.  **Chain Assignment**: Iterates through fragments to assign Chain IDs.
7.  **Merging**: Combines fragments into a single AtomGroup.
8.  **Transformation**: Unwraps, centers, and wraps coordinates to fix PBC.
9.  **Output**: Writes cleaned structure files (.pdb, .gro) and the transformed
    trajectory (.xtc).

Parameters (Global Variables)
-----------------------------
path : str
    Base directory for the project code.
step : int
    Stride used for slicing the output trajectory.
last_frames : int
    Number of frames to include from the end of the simulation.
output_path : str
    Directory containing the simulation outputs (must be set by user).

Outputs
-------
Files are written to the `postprocessing` subdirectory of `output_path`:
* `*_protein_pi.pdb`: Initial structure of protein and ions.
* `*_pi_clean.pdb`: Processed structure (PDB format).
* `*_pi_clean.gro`: Processed structure (GRO format).
* `*_full_pi_ref.xtc`: Transformed trajectory file.
"""

import sys
import os
import glob
import pickle
import argparse
from string import ascii_uppercase
import numpy as np
import MDAnalysis as mda
from MDAnalysis import transformations
from MDAnalysis.coordinates.PDB import PDBWriter


def get_args(argv=None):

    parser = argparse.ArgumentParser(description='input_output')
    #parser.add_argument('--path_tools', action="store", dest='path_tools', default='./')
    #parser.add_argument('--path', action="store", dest='path', default='./')
    parser.add_argument('--output', action="store", dest='output', default='./')
    parser.add_argument('--label', action="store", dest='label', default='run')
    #parser.add_argument('--result', action="store", dest='results_name', default='result_files')
    #parser.add_argument('--pwi', action="store", dest='bool_pwi', default='False')
    parser.add_argument('--step', action="store", dest='step', default=1)
    parser.add_argument('--frames', action="store", dest='frames', default=0)

    return parser.parse_args(argv)
 

def run_transform(output_path, label, step, last_frames): 

    # #### Load Metadata
    
    bookkeeping_path = os.path.join(output_path, 'postprocessing')
    print("bookkeeping_path", bookkeeping_path)
    pattern = '*_topol_mix.top'
    files = glob.glob(os.path.join(bookkeeping_path, pattern))
    label_new = files[0].split('/')[-1].strip().replace('_topol_mix.top', '')
    
    top = bookkeeping_path + "/" + label_new + "_md_run.tpr"
    
    
    #load all trajectory files
    with open(bookkeeping_path + '/' + label_new + '_xtc_list.pkl', 'rb') as file:
        xtcs = pickle.load(file)
    xtcs0 = xtcs[0]
    
    
    
    #select protein
    u1 = mda.Universe(top, xtcs0)
    prot1 = u1.select_atoms(
        "protein or resname ION or resname CL or resname NA "
        "and not name VWA VWB VWC VWD MW4"
    )
    with PDBWriter(bookkeeping_path + '/' + label_new + 'protein_pi.pdb') as pdb:
        pdb.write(prot1)
    
    
    # select fragments
    u = mda.Universe(top, xtcs, continuous=False, in_memory=False)
    prot = u.select_atoms(
        "protein or resname ION or resname CL or resname NA "
        "and not name VWA VWB VWC VWD MW4"
    )
    
    # add chain IDS inluding ions
    u.add_TopologyAttr("chainID")
    fragments = [
        f for f in prot.fragments
        if len(f) > 1 or set(f.atoms.residues.resnames).intersection(
            {'ION', 'CL', 'NA'}
        )
    ]
    chain_ids = ascii_uppercase  # A-Z
    for i, frag in enumerate(fragments):
        chain_id = chain_ids[i % len(chain_ids)]
        frag.atoms.chainIDs = chain_id
    
    # make one atom group out of the fragements in fragments
    all_frag_atoms = fragments[0]
    for frag in fragments[1:]:
        all_frag_atoms += frag
    
    all_frag_atoms.residues.resids = np.arange(
        1, len(all_frag_atoms.residues) + 1
    )
    
    # write including fragment info
    all_frag_atoms.unwrap(compound="fragments")
    all_frag_atoms.write(bookkeeping_path + '/' + label_new + "_pi_clean.pdb")
    all_frag_atoms.write(bookkeeping_path + '/' + label_new + '_pi_clean.gro')
    
    
    # do transforms
    try:
        workflow_ref = [
            transformations.unwrap(all_frag_atoms),
            transformations.center_in_box(all_frag_atoms.fragments[0]),
            transformations.wrap(all_frag_atoms, compound='fragments')
        ]
    
    except Exception:
        workflow_ref = [
            transformations.unwrap(all_frag_atoms),
            transformations.center_in_box(all_frag_atoms.fragments[0])
        ]
    
    all_frag_atoms.universe.trajectory.add_transformations(*workflow_ref)
    
    # write to file
    if last_frames==0:
        all_frag_atoms.write(
            bookkeeping_path + '/' + label_new + '_full_pi_ref.xtc',
            frames=all_frag_atoms.universe.trajectory[::step]
        )
    else:
        all_frag_atoms.write(
            bookkeeping_path + '/' + label_new + '_full_pi_ref.xtc',
            frames=all_frag_atoms.universe.trajectory[-last_frames::step]
        )

    return 0


if __name__ == "__main__":
    
    args = get_args()
    run_transform( args.output, args.label, int(args.step), int(args.frames))