import re
import os
import shutil
import glob
import itertools
import numpy as np
import MDAnalysis as mda
import numpy as np
import math
from scipy.spatial.distance import cdist


def calc_centroid(coords):
    """
    Calculate the geometric center (centroid) of a set of coordinates.

    Parameters
    ----------
    coords : np.ndarray
        Array of shape (N, 3) containing atomic coordinates.

    Returns
    -------
    np.ndarray
        Array of shape (3,) representing the centroid [x, y, z].
    """
    return np.mean(coords, axis=0)

def calc_normal(a, b, c):
    """
    Calculate the unit normal vector to a plane defined by three points.

    Parameters
    ----------
    a : np.ndarray
        Coordinates of point A (3,).
    b : np.ndarray
        Coordinates of point B (3,).
    c : np.ndarray
        Coordinates of point C (3,).

    Returns
    -------
    np.ndarray
        Normalized vector (3,) perpendicular to the plane formed by A, B, and C.
    """
    # Normal to plane defined by 3 points
    v1, v2 = b - a, c - a
    n = np.cross(v1, v2)
    return n / np.linalg.norm(n)

def calc_angle(v1,v2):
    """
    Calculate angles between a single vector v1 and a set of vectors v2.

    Parameters
    ----------
    v1 : np.ndarray
        Reference vector of shape (3,).
    v2 : np.ndarray
        Array of vectors of shape (N, 3).

    Returns
    -------
    np.ndarray
        Array of shape (N,) containing angles in degrees.
    """
    v1_u = v1 / np.linalg.norm(v1)
    # 3. Normalize v2 (Vectorized)
    # axis=1 computes norm per row. 
    # keepdims=True ensures shape is (N, 1) so we can divide (N, 3) by it.
    v2_norms = np.linalg.norm(v2, axis=1, keepdims=True)
    v2_u = v2 / v2_norms
    
    # 4. Dot Product (Matrix-Vector Multiplication)
    # v2_u is (N, 3), v1_u is (3,). Result is (N,)
    dots = np.dot(v2_u, v1_u)
    
    # 5. Clip and Calculate Angle
    return np.degrees(np.arccos(np.clip(dots, -1.0, 1.0)))


def calc_angle_scalar(v1, v2):
    """
    Calculate the angle between two individual vectors.

    Parameters
    ----------
    v1 : np.ndarray
        Vector of shape (3,).
    v2 : np.ndarray
        Vector of shape (3,).

    Returns
    -------
    float
        Angle between v1 and v2 in degrees.
    """
    v1u, v2u = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
    return np.degrees(np.arccos(np.clip(np.dot(v1u, v2u), -1.0, 1.0)))

def get_interaction_sel(key, group_a, group_b):
    """
    Identify potential interacting residue pairs based on interaction type.

    Parameters
    ----------
    key : str
        Type of interaction ('cation_pi', 'pi_stacking', 'hbond', 'salt_bridge').
    group_a : MDAnalysis.core.groups.AtomGroup
        First selection of atoms.
    group_b : MDAnalysis.core.groups.AtomGroup
        Second selection of atoms.

    Returns
    -------
    np.ndarray
        Array of unique residue index pairs (N, 2) relevant for the interaction.
    """
    l_sel=[]

    if key=='cation_pi':
                
        str1="(resname ARG and (name NH1 NH2)) or (resname LYS and name NZ)"
        str2="resname PHE TYR TRP and (name CG CE1 CE2 CD2 CZ2 CZ3)"
        
        valid_type_pairs= [('LYS', 'PHE'), ('ARG', 'PHE'), ('LYS', 'TYR'), ('ARG', 'TYR')] #cations, aromatics
        
        l_valid_pairs=[]
        
        for t_1, t_2 in valid_type_pairs:
            
            #cations
            sel_a1 = group_a.select_atoms(f"{str1} and resname {t_1}")
            sel_b1 = group_b.select_atoms(f"{str1} and resname {t_1}")
            
            #aromatics
            sel_a2 = group_a.select_atoms(f"{str2} and resname {t_2}")
            sel_b2 = group_b.select_atoms(f"{str2} and resname {t_2}")

            #combination (cation, aromatics)
            p1=itertools.product((sel_a1).resids, (sel_b2).resids)
            p2=itertools.product((sel_b1).resids, (sel_a2).resids)
        
            #save pairs
            l_valid_pairs.append(p1)
            l_valid_pairs.append(p2)
            
        final_valid_res_pairs = itertools.chain.from_iterable(l_valid_pairs)
        #print(type(final_valid_res_pairs))
        
        
        
    if key=='pi_stacking':
        
        str1="resname PHE TYR TRP and (name CG CE1 CE2 CD2 CZ2 CZ3)"

        sel_a = group_a.select_atoms(str1)
        sel_b = group_b.select_atoms(str1)

        final_valid_res_pairs = itertools.chain.from_iterable([itertools.product((sel_a).resids, (sel_b).resids)])


    
    if key=='hbond':

        str1="protein and (name N NE ND1 ND2 NE2 NZ OG OG1 OH)"
        str2="protein and (name O OD1 OD2 OE1 OE2 OG OG1 OH ND1 NE2)"

        valid_type_pairs= [('ARG', 'ASP'), ('ARG', 'GLU'), ('LYS', 'ASP'), ('LYS', 'GLU')] #to be done


        l_valid_pairs=[]
        
        for t_1, t_2 in valid_type_pairs:
            
            #donor
            sel_a1 = group_a.select_atoms(f"{str1} and resname {t_1}")
            sel_b1 = group_b.select_atoms(f"{str1} and resname {t_1}")
            

            #acceptor
            sel_a2 = group_a.select_atoms(f"{str2} and resname {t_2}")
            sel_b2 = group_b.select_atoms(f"{str2} and resname {t_2}")
 
        
            #combination (donor, acceptor)
            p1=itertools.product(sel_a1, sel_b2)
            p2=itertools.product(sel_b1, sel_a2)
        
            #save pairs
            l_valid_pairs.append(p1)
            l_valid_pairs.append(p2)
        
        final_valid_res_pairs = itertools.chain.from_iterable(l_valid_pairs)
       
        
    if key=='salt_bridge':


        str1="(resname ARG and (name NH1 NH2)) or (resname LYS and name NZ)"
        str2="(resname ASP and (name OD1 OD2)) or (resname GLU and (name OE1 OE2))"


        valid_type_pairs= [('ARG', 'ASP'), ('ARG', 'GLU'), ('LYS', 'ASP'), ('LYS', 'GLU')] #donor, acceptor

          
        l_valid_pairs=[]
        
        for t_1, t_2 in valid_type_pairs:
            
            #donor
            sel_a1 = group_a.select_atoms(f"{str1} and resname {t_1}")
            sel_b1 = group_b.select_atoms(f"{str1} and resname {t_1}")
            

            #acceptor
            sel_a2 = group_a.select_atoms(f"{str2} and resname {t_2}")
            sel_b2 = group_b.select_atoms(f"{str2} and resname {t_2}")
 
        
            #combination (donor, acceptor)
            p1=itertools.product((sel_a1).resids, (sel_b2).resids)
            p2=itertools.product((sel_b1).resids, (sel_a2).resids)
        
            #save pairs
            l_valid_pairs.append(p1)
            l_valid_pairs.append(p2)
        
        final_valid_res_pairs = itertools.chain.from_iterable(l_valid_pairs)
        
    return np.unique(np.array(list(final_valid_res_pairs)), axis=0)

# --- Main cation–π detection function ---
def cation_pi_contact(candidate_res_pairs, distance_cutoff=6.0, angle_cutoff=60.0):
    """
    Detect cation–π contacts between cations and aromatic residues.

    Parameters
    ----------
    candidate_res_pairs : iterator
        Iterable of (cations, aromatics) tuples where elements are MDAnalysis Residue objects.
    distance_cutoff : float, optional
        Max distance between cation atom and aromatic centroid (default 5.0).
    angle_cutoff : float, optional
        Max deviation from normal vector in degrees (default 30.0).

    Returns
    -------
    list of list
        List of [cation_resid, aromatic_resid] pairs satisfying criteria.
    """
    #print('aromatics_a', aromatics.atoms.names, len(aromatics.atoms), flush=True)
    #print('cations',cations.atoms.names, len(cations.atoms), flush=True)

    aromatic_atoms = {
    "PHE": ["CG", "CE1", "CE2"],
    "TYR": ["CG", "CE1", "CE2"],
    "TRP": ["CD2", "CZ2", "CZ3"],
    }
    l_res_pairs=[]
    for cations, aromatics in candidate_res_pairs: 
        names_array = aromatics.atoms.names

        # Create a boolean mask (True where name is in your list)
        mask = np.isin(names_array, aromatic_atoms[aromatics.resname])
        ring= aromatics.atoms[mask]

        if len(ring) < 3:
            continue
        
        coords =np.array(ring.positions)
        centroid = calc_centroid(coords)
        #print("centr", centroid)
        normal = calc_normal(*coords[:3])
        #print("normal", normal)
        vec = cations.atoms.positions - centroid #for all atoms in the cation
        mask_dist=np.linalg.norm(vec, axis=1)<= distance_cutoff
        
        if sum(mask_dist)>0:
            angles=calc_angle(normal, vec)
            #print("angles", angles)
            min_angles=np.minimum(abs(angles), abs(180 - angles))
            #print("min_angles", min_angles)
            mask_angle=min_angles<= angle_cutoff
            if sum(mask_angle)>0:
                l_res_pairs.append([int(cations.resid), int(aromatics.resid)])

    return [list(x) for x in set(tuple(p) for p in l_res_pairs)]


# --- Main π–stacking detector ---

def pi_stacking_contact(candidate_res_pairs,
                        distance_cutoff=7.0, angle_cutoff=30.0, psi_cutoff=45.0):
    """
    Detect π–stacking contacts between aromatic residues.

    Parameters
    ----------
    candidate_res_pairs : iterator
        Iterable of (aromatics_a, aromatics_b) tuples.
    distance_cutoff : float, optional
        Max distance between centroids (default 6.0).
    angle_cutoff : float, optional
        Max angle between plane normals (default 30.0).
    psi_cutoff : float, optional
        Min angle between normal and centroid vector (default 60.0).

    Returns
    -------
    list of list
        List of [resid_a, resid_b] pairs satisfying criteria.
    """    
    aromatic_atoms = {
    "PHE": ["CG", "CE1", "CE2"],
    "TYR": ["CG", "CE1", "CE2"],
    "TRP": ["CD2", "CZ2", "CZ3"],
}
    #print('aromatics_a', aromatics_a.atoms.names, flush=True)
    #print('aromatics_b',aromatics_b.atoms.names, flush=True)

    l_res_pairs=[]
    for aromatics_a, aromatics_b in candidate_res_pairs: 
        names_array_a = aromatics_a.atoms.names
        # Create a boolean mask (True where name is in your list)

        #a
        mask_a = np.isin(names_array_a, aromatic_atoms[aromatics_a.resname])
        ring_a= aromatics_a.atoms[mask_a]
        if len(ring_a) < 3:
            continue
            
        coords_a =np.array(ring_a.positions)
        centroid_a = calc_centroid(coords_a)       
        normal_a = calc_normal(*coords_a[:3])

        #b
        names_array_b = aromatics_b.atoms.names
        mask_b = np.isin(names_array_b, aromatic_atoms[aromatics_b.resname])
        ring_b= aromatics_b.atoms[mask_b]
        if len(ring_b) < 3:
            continue

        coords_b =np.array(ring_b.positions)
        centroid_b = calc_centroid(coords_b)
        normal_b = calc_normal(*coords_b[:3])

        #print('cent_a',centroid_a)
        #print('cent_b',centroid_b)
        dist = np.linalg.norm(centroid_b - centroid_a, axis=0)
        #print(dist)
        mask_dist=dist <= distance_cutoff
        #print('dist',dist)
        if mask_dist:
            #print(normal_a)
            #print(normal_b)
            angle_norm = calc_angle_scalar(normal_a,normal_b)
            #print('angle norm', angle_norm)
            angle=np.abs(np.minimum(angle_norm,180-angle_norm))
            mask_angle=angle <= angle_cutoff
            #print('angle', angle)

            if mask_angle:
                centroid_vec_a = centroid_b - centroid_a
                centroid_vec_b = centroid_a - centroid_b
                psi_a = calc_angle_scalar(normal_a, centroid_vec_a)
                psi_min_a = np.abs(np.minimum(psi_a, 180- psi_a ))
                psi_b = calc_angle_scalar(normal_b, centroid_vec_b)
                psi_min_b = np.abs(np.minimum(psi_b, 180 - psi_b ))
                psi_min = np.minimum(psi_min_a, psi_min_b)
                #print('psi_min', psi_min)
                if psi_min <= psi_cutoff:
                    l_res_pairs.append([int(aromatics_a.resid), int(aromatics_b.resid)])
            #return True

    return [list(x) for x in set(tuple(p) for p in l_res_pairs)]





def hbond_contact(donors, hydrogens, acceptors, distance_cutoff=3.5, angle_cutoff=30.0):
    """
    Detect hydrogen bonds between donors and acceptors using specific hydrogens.

    Parameters
    ----------
    donors : AtomGroup
        MDAnalysis AtomGroup containing donor atoms.
    hydrogens : AtomGroup
        MDAnalysis AtomGroup containing hydrogen atoms.
    acceptors : AtomGroup
        MDAnalysis AtomGroup containing acceptor atoms.
    distance_cutoff : float, optional
        Max distance between donor and acceptor (default 3.5).
    angle_cutoff : float, optional
        Max deviation from linear (180 deg) for D-H...A angle (default 30.0).

    Returns
    -------
    list of list
        List of [donor_resid, acceptor_resid] pairs.
    """

    #print('donors', donors.atoms.names, len(donors.atoms), flush=True)
    #print('hydrogens',hydrogens.atoms.names, len(hydrogens.atoms), flush=True)
    #print('acceptors',acceptors.atoms.names, len(acceptors.atoms), flush=True)

    
    l_res_pairs=[]
    for donor in donors:
        for hydrogen in hydrogens:
            # Only consider hydrogens near the donor (likely bonded)
            if np.linalg.norm(hydrogen.position - donor.position) > 1.2:
                continue

            for acceptor in acceptors:
                if acceptor.residue == donor.residue:
                    continue

                # Distance check (D···A)
                dist = np.linalg.norm(donor.position - acceptor.position)
                if dist > distance_cutoff:
                    continue

                # Angle check (D–H···A)
                v1 = hydrogen.position - donor.position
                v2 = acceptor.position - hydrogen.position
                angle = np.degrees(np.arccos(
                    np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
                ))

                if angle >= (180 - angle_cutoff):
                    l_res_pairs.append([int(donor.residue.resid), int(acceptor.residue.resid)])
            #return True

    #return [list(x) for x in set(tuple(p) for p in l_res_pairs)]
    return [list(x) for x in [tuple(p) for p in l_res_pairs]]


def salt_bridge_contact(candidate_res_pairs, distance_cutoff=4.0):
    """
    Detect salt bridges between cationic and anionic residues.

    Parameters
    ----------
    candidate_res_pairs : iterator
        Iterable of (cat_res, an_res) tuples where elements are Residue objects.
    distance_cutoff : float, optional
        Maximum heavy-atom distance (Å) to define a salt bridge (default 4.0).

    Returns
    -------
    list of list
        List of [cat_resid, an_resid] pairs.
    """

    #print('anions', anions.atoms.names, len(anions.atoms), flush=True)
    #print('cations',cations.atoms.names, len(cations.atoms), flush=True)

    
    l_res_pairs=[]

    #filtering for ab, ba pairs that are unpysical (will never make a salt bridge)
    for cat_res, an_res in candidate_res_pairs:

        #dist = np.linalg.norm(np.array(cat_res.atoms.positions) - np.array(an_res.atoms.positions), axis=1)
        dist_matrix = cdist(cat_res.atoms.positions, an_res.atoms.positions, metric='euclidean')

        #print(dist_matrix)

        exists = np.any(dist_matrix <= distance_cutoff)

        if exists:
            l_res_pairs.append([int(cat_res.resid), int(an_res.resid)])
            
    return [list(x) for x in set(tuple(p) for p in l_res_pairs)]
