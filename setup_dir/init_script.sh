#!/bin/bash

# exit when any command fails
set -e

# git_repo=$1;
git_repo="https://github.com/CarlosJesusGH/vina_docking.git"
DIR="/home/jovyan/vina_docking/"
if [ -d "$DIR" ]; then
### Take action if $DIR exists ###
cd $DIR
git reset --hard origin/master
git pull
else
###  Control will jump here if $DIR does NOT exists ###
git clone $git_repo
# mv repo_dir new_repo_dir
fi
# instead of conda init + logout/login
# source /root/miniconda3/etc/profile.d/conda.sh
# conda deactivate
# conda activate Env