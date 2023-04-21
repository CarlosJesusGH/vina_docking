# set docking parameters
receptor_molecule = 'protein.pdbqt'
vina_map_center = [-57.735, -3.936, -24.746]
vina_map_box_size = [41.25, 47.25, 45.0]
docking_exhaustiveness = 10 #32
docking_n_poses = 5 #20

# read files in directory according to pattern
# this could be changed if we know the exact name of the ligand files
import glob
ligand_molecules = sorted(glob.glob('drug*.pdbqt'))
print(ligand_molecules)

from vina import Vina
import time

v = Vina(sf_name='vina')

v.set_receptor(receptor_molecule)

# set a time flag for each ligand
global_init_time = time.time()

for i, ligand_molecule in enumerate(ligand_molecules):
  # logging
  print("working with ligand %d out of %d" % (i, len(ligand_molecules)))
  print("name of ligand file is %s" % ligand_molecule)
  ligand_init_time = time.time()

  v.set_ligand_from_file(ligand_molecule)
  v.compute_vina_maps(center=vina_map_center, box_size=vina_map_box_size)

  # Score the current pose
  energy = v.score()
  print('Score before minimization: %.3f (kcal/mol)' % energy[0])

  # Minimized locally the current pose
  energy_minimized = v.optimize()
  print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
  v.write_pose(ligand_molecule.split('.')[0] + '_minimized.pdbqt', overwrite=True)

  # Dock the ligand
  v.dock(exhaustiveness=docking_exhaustiveness, n_poses=docking_n_poses)
  v.write_poses(ligand_molecule.split('.')[0] + '_vina_out.pdbqt', n_poses=5, overwrite=True)

  # log time in human readable format
  print("time elapsed for ligand is %s" % time.strftime("%H:%M:%S", time.gmtime(time.time() - ligand_init_time)))
  print("global time elapsed is %s" % time.strftime("%H:%M:%S", time.gmtime(time.time() - global_init_time)))
  print("-"*50)