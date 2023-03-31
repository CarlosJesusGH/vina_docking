FROM jupyter/minimal-notebook:python-3.7
MAINTAINER Carlos Garcia-Hernandez carlos.garcia2@bsc.es

# -----------------------------------------------------------------------
# VARIABLE DECLARATION

ENV HOME_DIR=/home
ENV REPO_DIR=${HOMEDIR}/vina_docking
ENV SETUP_DIR=${REPODIR}/setup_dir
ENV SERVER_DIR=${REPODIR}/tools
ENV INIT_DIR=${REPODIR}/init_dir

# -----------------------------------------------------------------------
# UPDATE SYSTEM AND CLONE GIT REPOSITORY

RUN apt -y update && apt -y upgrade
RUN apt install -y git
WORKDIR ${HOME_DIR}
RUN git clone https://github.com/CarlosJesusGH/vina_docking.git
RUN ls ${REPO_DIR}

# -----------------------------------------------------------------------
# SETUP THE SYSTEM

RUN apt install -y wget nano tree htop ncdu curl

WORKDIR ${SETUP_DIR}


RUN echo "-----------------------------------------------------------------------"

# Initialize conda in bash config fiiles:
RUN conda init
RUN conda init bash

RUN conda --version
RUN conda env list

RUN echo "-----------------------------------------------------------------------"

# Create GC3Env conda environment and install requirements
RUN conda create -y -n GC3Env python=2.7

# Make RUN commands use the new environment:    [1]
# SHELL ["conda", "run", "-n", "base", "/bin/bash", "-c"]

# Make RUN commands use the new environment:    [1]
SHELL ["conda", "run", "-n", "GC3Env", "/bin/bash", "-c"]
RUN echo $CONDA_DEFAULT_ENV
RUN python --version
RUN \
    pip install -r requirements_GC3Env.txt && \
    # pip install mysql-python
    pip freeze

RUN echo "-----------------------------------------------------------------------"

RUN \
    pip install MySQL-python==1.2.5 && \
    pip install scikit-image==0.11.3 && \
    pip freeze

RUN echo "-----------------------------------------------------------------------"

# to start the service, maybe we could also use something like this: [3]
# SHELL ["/sbin/init", "-d"]
# I didn't try it but it might work

# The default shell on Linux is ["/bin/sh", "-c"] [2]
SHELL ["/bin/sh", "-c"]

RUN apt install -y mysql-server
RUN apt install -y apache2 apache2-utils

RUN ls /var/lib/mysql/
RUN ls /etc/mysql/

RUN service mysql start
RUN service apache2 start

RUN \
    service mysql restart && \
    mysql --execute="CREATE DATABASE graphcrunch3; CREATE DATABASE gc3;" && \
    mysql --execute="SHOW DATABASES;"

RUN echo "-----------------------------------------------------------------------"

# RUN service mysql restart && mysql --execute="SHOW DATABASES;"

# RUN tar -xzf WebServer_bkp_*.tar.gz -C ${SERVERDIR}
RUN pwd && ls
WORKDIR ${SERVERDIR}

SHELL ["conda", "run", "-n", "GC3Env", "/bin/bash", "-c"]
RUN \
    service mysql restart && \
    ./manage.py makemigrations --setting=WebServer.settings_prod && \
    ./manage.py migrate --setting=WebServer.settings_prod --noinput && \
    ./manage.py syncdb --setting=WebServer.settings_prod --noinput && \
    python manage.py migrate

RUN echo "-----------------------------------------------------------------------"
# <!-- create environment with python 3.7 and install all the necessary packages -->

WORKDIR ${SETUPDIR}

SHELL ["/bin/sh", "-c"]
# RUN conda deactivate
RUN conda update conda -y
RUN conda config --append channels conda-forge
RUN conda config --append channels bioconda
#conda create --name env_37 --file /home/Downloads/GC3-WWW/drive_files/requirements_env37.txt
RUN conda create --name env_37 python=3.7
# RUN conda activate env_37
# SHELL ["conda", "run", "-n", "env_37", "/bin/bash", "-c"]
# RUN pwd && ls
RUN conda install --name env_37 --yes --file requirements_env37v2.txt
#conda install -c intel -y mkl-fft==1.3.0 mkl-random==1.2.2
# RUN conda install --name env_37 --yes -c intel -y mkl-fft mkl-random
RUN conda install -y --name env_37 -c conda-forge mkl_fft mkl_random
RUN conda install --name env_37 --yes -c anaconda setuptools ipython_genutils
RUN conda install --name env_37 --yes pandas==1.3.4

RUN echo "-----------------------------------------------------------------------"
# <!-- create environment with python 2.7 and install all the necessary packages -->

# conda deactivate
#conda create --name env_27 python=2.7 numpy networkx matplotlib ipykernel scipy scikit-learn -y
RUN conda create --name env_27 python=2.7 numpy networkx==1.11 matplotlib ipykernel scipy scikit-learn -y
RUN conda install -n env_27 -y -c bioconda gnuplot
# cp /home/Downloads/GC3-WWW/drive_files/icell_input_files/gnuplot-py-1.8.tar.gz .
RUN tar -xf gnuplot-py-1.8.tar.gz
# RUN rm gnuplot-py-1.8.tar.gz
RUN conda install -n env_27 decorator=4.4.0 -y
# conda deactivate
SHELL ["conda", "run", "-n", "env_27", "/bin/bash", "-c"]
# source activate env_27 && 
RUN python --version && cd gnuplot-py-1.8 && python setup.py install
# source activate env_27 && 
# RUN python --version
SHELL ["conda", "run", "-n", "env_27", "python", "-c"]
RUN import networkx; import Gnuplot; print('hello gnuplot world');
# cd ..
SHELL ["/bin/sh", "-c"]
RUN pwd && ls
RUN rm -r ${SETUPDIR}

RUN echo "-----------------------------------------------------------------------"

WORKDIR ${INITDIR}

EXPOSE 8000

RUN cp ${STARTDIR}/init_script.sh .
RUN chmod 744 init_script.sh

RUN rm -r ${REPODIR}


ENTRYPOINT ["/home/init_script.sh"]
