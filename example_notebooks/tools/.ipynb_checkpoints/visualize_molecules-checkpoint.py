import py3Dmol
from rdkit import Chem

def plot_3d_ligand(mol):
    # mol_raw = Chem.SDMolSupplier(ligand_path)[0]
    # mol_prep = ligand_preparation.prepare_mols_from_sdf(ligand_path)[0]
    # --------------
    
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    
    view.addModel(Chem.MolToMolBlock(mol,False),'mol')
    x = view.getModel()
    x.setStyle({},{'stick':{'colorscheme':'cyanCarbon','radius':0.2}})
    
    view.zoomTo()
    view.show()
    
def plot_prot_and_ligand(prot_path, lig_str, include_box=False, center=[], size=[], add_surface=False, mol=None):
    view = py3Dmol.view()
    view.removeAllModels()
    view.setViewStyle({'style':'outline','color':'black','width':0.1})
    
    # add receptor
    view.addModel(open(prot_path,'r').read(),format='pdb')
    Prot=view.getModel()
    Prot.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    if add_surface:
        view.addSurface(py3Dmol.VDW,{'opacity':0.6,'color':'white'})
    
    # add ligand
    # view.addModel(lig_str,format='pdb')
    # ref_m = view.getModel()
    # ref_m.setStyle({'cartoon':{'arrows':True, 'tubes':True, 'style':'oval', 'color':'white'}})
    # ref_m.setStyle({'stick':{'style':'oval', 'color':'white'}})
    # ref_m.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
    
    # print("lig_str", lig_str)
    # results=Chem.SDMolSupplier(lig_str)
    # p=Chem.MolToMolBlock(results[0],False)
    
    view.addModel(Chem.MolToMolBlock(lig_str,False),'mol')
    ligModel = view.getModel()
    ligModel.setStyle({},{'stick':{'colorscheme':'magentaCarbon','radius':0.2}})
    
    if include_box:
        view.addBox({"center":{"x":center[0],"y":center[1],"z":center[2]},"dimensions":{"w":size[0],"h":size[1],"d":size[2]},"color":'red','opacity':0.3});
    
    view.zoomTo()
    view.show()