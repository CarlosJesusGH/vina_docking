import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from scipy.spatial import distance
from meeko import MoleculePreparation

from . import parse_molecules

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
    
def translate_ligand_from_meekoprep(meeko_prep, new_center):
    # translate ligand to new center
    for atom_index in meeko_prep.setup.coord:
        x, y, z = meeko_prep.setup.coord[atom_index]
        meeko_prep.setup.coord[atom_index] = np.array([x + new_center[0], y + new_center[1], z + new_center[2]])
   
# TODO
# def translate_ligand_from_mol(mol, new_center):
#     pass

def translate_ligand_from_pdbqt(pdbqt_path, new_center, verbose=True):
    # lig = parse_molecules.MolFromPDBQTBlock(block=open(pdbqt_path,'r').read(),sanitize=True,removeHs=False)
    lig = parse_molecules.MolFromPDBQTBlock(block=pdbqt_path,sanitize=True,removeHs=False)
    protonated_lig = rdkit.Chem.AddHs(lig)
    # generate 3D coordinates for the ligand
    rdkit.Chem.AllChem.EmbedMolecule(protonated_lig)
    # convert to a PDBQT string using the MoleculePreparation class from Meeko.
    meeko_prep = MoleculePreparation()
    meeko_prep.prepare(protonated_lig)
    # preparator.show_setup()
    if verbose: print("initial_ligand_center", get_ligand_center(meeko_prep))
    translate_ligand_from_meekoprep(meeko_prep, new_center)
    if verbose: print("translated_ligand_center", get_ligand_center(meeko_prep))
    # At this point, pdbqt_string can be written to a file for docking with AutoDock-GPU or Vina, or passed directly to Vina within Python using set_ligand_from_string(pdbqt_string)
    # with open("./aux_pdbqt.pdbqt", "w") as text_file:
    #     text_file.write(mkprep.write_pdbqt_string())
    lig_pdbqt_string = meeko_prep.write_pdbqt_string()
    return lig_pdbqt_string
        
def get_optimal_box_size_for_ligand(ligand_path):
    import subprocess as sp
    output = sp.getoutput('perl /home/jovyan/vina_docking/tools/3dparty_scripts/eBoxSize-1.1.pl ' + ligand_path)
    return output
    