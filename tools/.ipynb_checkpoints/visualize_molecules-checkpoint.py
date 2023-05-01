import py3Dmol
from rdkit import Chem
import ast
   
def show_3d_ligand(mol, size=(500, 400), style="stick", surface=False, opacity=0.5):
    """Draw molecule in 3D
    Args:
    ----
        mol: rdMol, molecule to show
        size: tuple(int, int), canvas size
        style: str, type of drawing molecule
               style can be 'line', 'stick', 'sphere', 'carton'
        surface, bool, display SAS
        opacity, float, opacity of surface, range 0.0-1.0
    Return:
    ----
        viewer: py3Dmol.view, a class for constructing embedded 3Dmol.js views in ipython notebooks.
    """
    assert style in ('line', 'stick', 'sphere', 'carton')
    mblock = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=size[0], height=size[1])
    viewer.addModel(mblock, 'mol')
    viewer.setStyle({style:{}})
    if surface:
        viewer.addSurface(py3Dmol.SAS, {'opacity': opacity})
    viewer.zoomTo()
    return viewer

def show_3d_receptor(receptor_path, receptor_name, size=(800, 600), 
                     show_pockets=False, pockets_opacity=0.5, 
                     show_drug_score=False,
                     show_surface=False, surface_opacity=0.5,
                    pocket_paths=None, pockets_data=None):
    view = py3Dmol.view(width=size[0], height=size[1])
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    view.addModel(open(receptor_path,'r').read(),'pdb')
    Prot=view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    if show_surface:
        view.addSurface(py3Dmol.VDW,{'opacity':surface_opacity,'color':'white'})
    if show_pockets:
        for cav_id, pocket_path in enumerate(pocket_paths, start=1):    
            drug_score = pockets_data["drug_score"].loc[cav_id]
            # Color from red to yellow according to drug_score (red-ish is better)
            color = ["#%02x%02x%02x" % (255, int(255*(1-drug_score)), 0)]
            view.addModel(open(pocket_path, 'r').read(), 'pqr')
            x = view.getModel()
            x.setStyle({},{'sphere':{'color':color[0],'opacity':pockets_opacity}})
            if show_drug_score:
                view.addLabel(str(drug_score), {'fontColor': color[0], 'position': {'x': 0, 'y': 0, 'z': 0}, 'backgroundColor': 'black', 'opacity': 0.7}, {'model': cav_id})
    view.zoomTo()
    return view
    
    
def show_3d_vina_score_dock(receptor_path, receptor_name, show_pockets, show_drug_score, show_surface, 
                            show_vina_box, ligand_style, rdkitmol_list, pocket_paths, pockets_data, scoring_row, 
                            show_only_best_pose):
    # Add the receptor
    view = show_3d_receptor(receptor_path, receptor_name, 
                            show_pockets=show_pockets, show_drug_score=show_drug_score, show_surface=show_surface, 
                            pocket_paths=pocket_paths, pockets_data=pockets_data)
    
    # Add the box
    if show_vina_box:
        center = ast.literal_eval(scoring_row['pocket_box_center'])
        size = ast.literal_eval(scoring_row['vina_box_size'])
        view.addBox({"center":{"x":center[0],"y":center[1],"z":center[2]},"dimensions":{"w":size[0],"h":size[1],"d":size[2]},"color":'red','opacity':0.3});
    
    # Add the ligands
    for i, mol in enumerate(rdkitmol_list):
        p=Chem.MolToMolBlock(mol,False)
        # print ('Pose: {} | Score: {}'.format(mol.GetProp('Pose'),mol.GetProp('Score')))
        view.addModel(p,'mol')
        x = view.getModel()
        if show_only_best_pose:
            x.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
            break
        # Color from red to yellow according to index
        color = "#%02x%02x%02x" % (255, int(255*(i/len(rdkitmol_list))), 0)
        x.setStyle({},{'stick':{'color':color,'radius':0.2}})
        
    # Center the view towards the ligand
    # view.zoomTo({'model':-1})
    # Re-center the viewer around the provided selection (unlike zoomTo, does not zoom).
    # ref: https://3dmol.csb.pitt.edu/doc/$3Dmol.GLViewer.html
    view.center({'model':-1})
    return view
    
    
# def smi2conf(smiles):
#     '''Convert SMILES to rdkit.Mol with 3D coordinates'''
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is not None:
#         mol = Chem.AddHs(mol)
#         AllChem.EmbedMolecule(mol)
#         AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
#         return mol
#     else:
#         return None