FROM jupyter/minimal-notebook:python-3.7
MAINTAINER Carlos Garcia-Hernandez carlos.garcia2@bsc.es

# -----------------------------------------------------------------------
# VARIABLE DECLARATION

ENV HOME_DIR=/home
ENV JOVYAN_DIR=${HOME_DIR}/jovyan
ENV REPO_DIR=${JOVYAN_DIR}/vina_docking
ENV SETUP_DIR=${REPO_DIR}/setup_dir

# -----------------------------------------------------------------------
# UPDATE SYSTEM AND CLONE GIT REPOSITORY

USER root

RUN apt -yq update && apt -yq upgrade
RUN apt install -yq git
WORKDIR ${HOME_DIR}
# TODO: change ADD for CLONE in the end
ADD . ${JOVYAN_DIR}
# RUN git clone https://github.com/CarlosJesusGH/vina_docking.git
RUN ls ${JOVYAN_DIR}


# -----------------------------------------------------------------------
# INSTALL OTHER USEFUL TOOLS

RUN apt install -yq wget nano tree htop ncdu curl
# RUN wget -O ${JOVYAN_DIR}/bsc_autodock_+_cli_+_python.ipynb 'https://docs.google.com/uc?export=download&id=1E19Clw-jJ3XtfLns9RINywS7qtZptWRF'
RUN \
    mv /usr/local/bin/start-notebook.sh /usr/local/bin/start-notebook_bkp.sh && \
    cp ${SETUP_DIR}/start-notebook.sh /usr/local/bin && \
    chmod 777 /usr/local/bin/start-notebook.sh

# necessary for schrodinger pymol-bundle
RUN apt install -yq libgl1

# -----------------------------------------------------------------------
# INSTALL ADFRsuite

WORKDIR ${SETUP_DIR}

RUN \
    wget https://ccsb.scripps.edu/adfr/download/1038/ -O ADFRsuite_x86_64Linux_1.0.tar.gz && \
    tar zxf ADFRsuite_x86_64Linux_1.0.tar.gz

RUN \
    replace_string="ans = ''\n    try:\n      ans = raw_input(text)\n    except EOFError as e:\n      print(e)" && \
    sed -i "s|ans = raw_input(text)|$replace_string|g" ./ADFRsuite_x86_64Linux_1.0/Tools/install.py

RUN \
    mkdir /ADFRsuite && \
    cd ADFRsuite_x86_64Linux_1.0/ && \
    ./install.sh -d /ADFRsuite -c 0

RUN \
    mv /ADFRsuite/bin/prepare_ligand /ADFRsuite/bin/prepare_ligand_bkp && \
    cp ${SETUP_DIR}/prepare_ligand /ADFRsuite/bin/ && \
    chmod 777 /ADFRsuite/bin/prepare_ligand

RUN echo "export PATH=/ADFRsuite/bin:\$PATH" >> ~/.bashrc
RUN echo "alias python=python3" >> ~/.bashrc
RUN echo "alias ipython=ipython3" >> ~/.bashrc

# -----------------------------------------------------------------------
# SHOW CONDA INFO

RUN conda --version && conda env list
# RUN source /opt/conda/etc/profile.d/conda.sh && conda activate base && conda env list

# -----------------------------------------------------------------------
# INSTALL VINA ON BASE CONDA ENVIRONMENT

SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]
RUN pip install vina
RUN pip install pdb-tools
RUN conda env list
# RUN conda list

# -----------------------------------------------------------------------
# SETUP MEEKO CONDA ENVIRONMENT (try doing it all on the same 'base' environment)

# RUN conda create -y --name env_meeko python=3.9
# SHELL ["conda", "run", "-n", "env_meeko", "/bin/bash", "-c"]
RUN conda install -c conda-forge numpy openbabel scipy rdkit -y
RUN pip install meeko
# RUN conda list

# -----------------------------------------------------------------------
# TOOLS TO VISUALIZE MOLECULES
RUN conda install -y -c conda-forge -c schrodinger pymol-bundle
RUN conda install -y -c conda-forge py3dmol mdanalysis fpocket
# prolif
 
# -----------------------------------------------------------------------
# DELETE UNNECESSARY FILES AND DIRS
RUN rm -rf ${REPO_DIR}/setup_dir
RUN rm -f ${REPO_DIR}/DocerfileAux.docker

# -----------------------------------------------------------------------
# STARTUP CONFIG - from: https://hub.docker.com/r/jupyter/base-notebook/dockerfile

EXPOSE 8888

# Configure container startup
ENTRYPOINT ["tini", "-g", "--"]
CMD ["start-notebook.sh", "--NotebookApp.token=''", "--NotebookApp.password=''", "--NotebookApp.allow_origin='*'"]

# SET EVERYTHING AS DEFAULT IMAGE
# WORKDIR ${JOVYAN_DIR}
WORKDIR $HOME
# Fix DL4006
# SHELL ["/bin/bash", "-o", "pipefail", "-c"]
# Fix permissions on /etc/jupyter as root
# USER root
# RUN fix-permissions /etc/jupyter/
# Switch back to jovyan to avoid accidental container runs as root
USER $NB_UID

# after any change, re-build and push new image to dockerhub:
# dockerbuild && dockerpush