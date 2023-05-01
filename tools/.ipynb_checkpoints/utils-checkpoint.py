# this script is from: https://github.com/AngelRuizMoreno/Jupyter_Dock/blob/main/utilities/utils.py

import py3Dmol
from pymol import cmd
from openbabel import pybel

from rdkit import Chem
from rdkit.Chem import AllChem,rdFMCS, rdMolAlign

# from pdbfixer import PDBFixer
# from openmm.app import PDBFile

import MDAnalysis as mda
from MDAnalysis.coordinates import PDB

import random, math

import numpy as np

def get_scores_from_pdbqt(pdbqt_file):
    results = [m for m in pybel.readfile(filename=pdbqt_file,format='pdbqt')]
    scores = []
    for pose in results:
        scores.append(float(pose.data['REMARK'].split()[2]))
    return scores


def getbox(selection='sele', extending = 6.0, software='vina'):
    
    ([minX, minY, minZ],[maxX, maxY, maxZ]) = cmd.get_extent(selection)

    minX = minX - float(extending)
    minY = minY - float(extending)
    minZ = minZ - float(extending)
    maxX = maxX + float(extending)
    maxY = maxY + float(extending)
    maxZ = maxZ + float(extending)
    
    SizeX = maxX - minX
    SizeY = maxY - minY
    SizeZ = maxZ - minZ
    CenterX =  (maxX + minX)/2
    CenterY =  (maxY + minY)/2
    CenterZ =  (maxZ + minZ)/2
    
    cmd.delete('all')
    
    if software == 'vina':
        return {'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}
    elif software == 'ledock':
        return {'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ}
    elif software == 'both':
        return ({'center_x':CenterX,'center_y': CenterY, 'center_z': CenterZ},{'size_x':SizeX,'size_y': SizeY,'size_z': SizeZ}),({'minX':minX, 'maxX': maxX},{'minY':minY, 'maxY':maxY}, {'minZ':minZ,'maxZ':maxZ})
    
    else:
        print('software options must be "vina", "ledock" or "both"')


def pdbqt_to_sdf(pdbqt_file=None,output=None, include_data=True):
    results = [m for m in pybel.readfile(filename=pdbqt_file,format='pdbqt')]
    out=pybel.Outputfile(filename=output,format='sdf',overwrite=True)
    for pose in results:
        if include_data:
            pose.data.update({'Pose':pose.data['MODEL']})
            pose.data.update({'Score':pose.data['REMARK'].split()[2]})
            del pose.data['MODEL'], pose.data['REMARK'], pose.data['TORSDO']
        out.write(pose)
    out.close()