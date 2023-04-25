import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from scipy.spatial import distance

def prepare_mols_from_sdf(filepath):
    mols = []
    for mol in Chem.SDMolSupplier(filepath):
        # check the protonation state of the molecule
        mol = Chem.AddHs(mol)
        # obtain coordinates for the molecule
        etkdg_params = rdDistGeom.ETKDGv3()
        rdDistGeom.EmbedMolecule(mol, etkdg_params)
        mols.append(mol)
    return mols

def get_ligand_center(meeko_prep):
    meeko_prep = meeko_prep.setup
    lig_xyz = []
    for atom_index, is_atom_ignored in meeko_prep.atom_ignore.items():
        if not is_atom_ignored:
            lig_xyz.append(meeko_prep.coord[atom_index].copy())
    lig_xyz = np.array(lig_xyz)
    # print("ligand center: %8.3f %8.3f %8.3f" % tuple(np.mean(lig_xyz, 0)))
    center = list(np.mean(lig_xyz, 0))
    # euclid_dist_origin = math.dist(center, [0,0,0])  # this works on python3.8+
    euclid_dist_origin = distance.euclidean(center, [0,0,0])
    return [ '%.2f' % elem for elem in center], euclid_dist_origin
    
def translate_ligand(meeko_prep, new_center):
    # translate ligand to new center
    for atom_index in meeko_prep.setup.coord:
        x, y, z = meeko_prep.setup.coord[atom_index]
        meeko_prep.setup.coord[atom_index] = np.array([x + new_center[0], y + new_center[1], z + new_center[2]])
        
def get_optimal_box_size_for_ligand(ligand_path):
    import subprocess as sp
    output = sp.getoutput('perl /home/jovyan/example_notebooks/tools/3dparty_scripts/eBoxSize-1.1.pl ' + ligand_path)
    return output
    