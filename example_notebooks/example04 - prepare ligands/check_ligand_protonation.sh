#!/bin/bash

source /opt/conda/etc/profile.d/conda.sh
conda activate env_meeko
# mk_prepare_ligand.py $@
python check_ligand_protonation.py $@
