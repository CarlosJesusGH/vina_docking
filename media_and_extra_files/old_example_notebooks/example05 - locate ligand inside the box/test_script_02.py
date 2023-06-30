
from vina import Vina
import numpy as np
from meeko import MoleculePreparation
from rdkit import Chem
from rdkit.Chem import rdDistGeom

def print_ligand_center(molsetup):
    lig_xyz = []
    for atom_index, is_atom_ignored in molsetup.atom_ignore.items():
        if not is_atom_ignored:
            lig_xyz.append(molsetup.coord[atom_index].copy())
    lig_xyz = np.array(lig_xyz)
    print("ligand center: %8.3f %8.3f %8.3f" % tuple(np.mean(lig_xyz, 0)))
    
def translate_ligand(mkprep, new_center):
    # translate ligand to new center
    for atom_index in mkprep.setup.coord:
        x, y, z = mkprep.setup.coord[atom_index]
        mkprep.setup.coord[atom_index] = np.array([x + new_center[0], y + new_center[1], z + new_center[2]])

# class vina.vina.Vina(sf_name='vina', cpu=0, seed=0, no_refine=False, verbosity=1)
v = Vina(sf_name='vina', seed=1)

v.set_receptor('1iep_receptor.pdbqt')

# pyridine = Chem.MolFromSmiles("C1=CC=CN=C1")
# pyridine = Chem.AddHs(pyridine)
# etkdg_params = rdDistGeom.ETKDGv3()
# rdDistGeom.EmbedMolecule(pyridine, etkdg_params)

mol = Chem.SDMolSupplier('1iep_ligand.sdf')[0]
mol = Chem.AddHs(mol)
etkdg_params = rdDistGeom.ETKDGv3()
rdDistGeom.EmbedMolecule(mol, etkdg_params)

mkprep = MoleculePreparation()
mkprep.prepare(mol)
# preparator.show_setup()
lig_string = mkprep.write_pdbqt_string()

v.set_ligand_from_string(lig_string)
maps_center=[15.190, 53.903, 16.917]
maps_size_angstroms=[20, 20, 20]
v.compute_vina_maps(center=maps_center, box_size=maps_size_angstroms)

# ---------

print_ligand_center(mkprep.setup)
translate_ligand(mkprep, maps_center)
print_ligand_center(mkprep.setup)

lig_string = mkprep.write_pdbqt_string()
v.set_ligand_from_string(lig_string)

score = v.score()
print("score:", score)
