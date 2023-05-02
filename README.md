# vina_docking
Tools to ease the automation of molecular docking using autodock vina and other related packages.

After installing Docker, execute this command to run our image:
```bash
docker run --rm -p 8888:8888 carlosjesusgh/vina_docking:latest
```

To create/save your own experiments, share a directory from your local machine to the Docker container
```bash
mkdir USER_EXPERIMENTS
docker run --rm -p 8888:8888 -v ./USER_EXPERIMENTS/:/home/jovyan/USER_EXPERIMENTS/ carlosjesusgh/vina_docking:latest
```

https://user-images.githubusercontent.com/8160204/235643111-be637b22-a382-47ea-ba35-47428eaba39e.mp4

Interactive 3D-view example

https://user-images.githubusercontent.com/8160204/235505816-cec032dd-5068-4e20-bc89-abe835047b67.mp4

To run your own example, open the notebook in 'test_experiments/real_case_step_by_step/', change the list of ligands and receptors according to your needs and click on 'Run/Run All Cells'

https://user-images.githubusercontent.com/8160204/235505977-8609263f-2be8-43bd-9582-4e5c335a43af.mp4
