from rdkit import Chem
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate

def MolFromPDBQTBlock(block, sanitize=True, removeHs=True):
    """Read PDBQT block to a RDKit Molecule
    Parameters
    ----------
        block: string
            Residue name which explicitly pint to a ligand(s).
        sanitize: bool (default=True)
            Should the sanitization be performed
        removeHs: bool (default=True)
            Should hydrogens be removed when reading molecule.
    Returns
    -------
        mol: rdkit.Chem.rdchem.Mol
            Molecule read from PDBQT
    """
    # mol = Chem.MolFromPDBBlock('\n'.join(pdb_lines), sanitize=False)
    pdbqt_mol = PDBQTMolecule.from_file(block, skip_typing=True)
    rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
    mol = rdkitmol_list[0]
    if sanitize:
        Chem.SanitizeMol(mol)
    else:
        Chem.GetSSSR(mol)
#     # reorder atoms using serial
#     new_order = sorted(range(mol.GetNumAtoms()),
#                        key=lambda i: (mol.GetAtomWithIdx(i)
#                                       .GetPDBResidueInfo()
#                                       .GetSerialNumber()))
#     mol = Chem.RenumberAtoms(mol, new_order)

#     # properties must be set on final copy of Mol, RenumberAtoms purges data
#     mol.SetProp('_Name', name)
#     for k, v in data.items():
#         mol.SetProp(str(k), str(v))

    return mol

# from: https://github.com/oddt/oddt
# file: https://github.com/oddt/oddt/blob/a3ff8b84b3abf986ad5bdbfebc9ef96cb8a84d8c/oddt/toolkits/extras/rdkit/__init__.py#L353
def MolFromPDBQTBlock_OLD(block, sanitize=True, removeHs=True):
    # this implementation is having problems with ATOM G used in vina
    """Read PDBQT block to a RDKit Molecule
    Parameters
    ----------
        block: string
            Residue name which explicitly pint to a ligand(s).
        sanitize: bool (default=True)
            Should the sanitization be performed
        removeHs: bool (default=True)
            Should hydrogens be removed when reading molecule.
    Returns
    -------
        mol: rdkit.Chem.rdchem.Mol
            Molecule read from PDBQT
    """
    pdb_lines = []
    name = ''
    data = {}
    for line in block.split('\n'):
        # Get all know data from REMARK section
        if line[:12] == 'REMARK  Name':
            name = line[15:].strip()
        elif line[:18] == 'REMARK VINA RESULT':
            tmp = line[19:].split()
            data['vina_affinity'] = tmp[0]
            data['vina_rmsd_lb'] = tmp[1]
            data['vina_rmsd_ub'] = tmp[2]

        # no more data to collect
        if line[:4] != 'ATOM':
            continue

        pdb_line = line[:56]
        pdb_line += '1.00  0.00           '

        # Do proper atom type lookup
        atom_type = line[71:].split()[1]
        if atom_type == 'A':
            atom_type = 'C'
        elif atom_type[:1] == 'O':
            atom_type = 'O'
        elif atom_type[:1] == 'H':
            atom_type = 'H'
            if removeHs:
                continue
        elif atom_type == 'NA':
            atom_type = 'N'

        pdb_lines.append(pdb_line + atom_type)
    mol = Chem.MolFromPDBBlock('\n'.join(pdb_lines), sanitize=False)
    if sanitize:
        Chem.SanitizeMol(mol)
    else:
        Chem.GetSSSR(mol)
    # reorder atoms using serial
    new_order = sorted(range(mol.GetNumAtoms()),
                       key=lambda i: (mol.GetAtomWithIdx(i)
                                      .GetPDBResidueInfo()
                                      .GetSerialNumber()))
    mol = Chem.RenumberAtoms(mol, new_order)

    # properties must be set on final copy of Mol, RenumberAtoms purges data
    mol.SetProp('_Name', name)
    for k, v in data.items():
        mol.SetProp(str(k), str(v))

    return mol