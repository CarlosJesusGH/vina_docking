# execute minimal-notebook with root user
docker run --rm -v /home/bscuser/vina_experiments:/home/jovyan/vina_experiments -p 8888:8888 --name vina_jupyter -e GRANT_SUDO=yes --user root jupyter/minimal-notebook:python-3.7 start-notebook.sh --NotebookApp.token='' --NotebookApp.password="" --NotebookApp.allow_origin="*"

# execute vina:no_tag
docker run --rm -v /home/bscuser/vina_experiments:/home/jovyan/vina_experiments -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root vina:no_tag start-notebook.sh --NotebookApp.token='' --NotebookApp.password='' --NotebookApp.allow_origin='*'

docker run -v /home/bscuser/vina_experiments:/home/jovyan/vina_experiments -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root vina:no_tag start-notebook-vina.sh --NotebookApp.token='' --NotebookApp.password='' --NotebookApp.allow_origin='*'

docker run --rm -v /home/bscuser/vina_experiments:/home/jovyan/vina_experiments -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root vina:no_tag "source ~/.bashrc && start-notebook.sh --NotebookApp.token='' --NotebookApp.password='' --NotebookApp.allow_origin='*'"


docker run --rm -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root vina:no_tag start-notebook.sh --NotebookApp.token='' --NotebookApp.password='' --NotebookApp.allow_origin='*'

docker run --rm -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root vina:no_tag

# execute carlosjesusgh/vina_docking:v23MMDD
docker run --rm -p 8888:8888 --name vina_cont -e GRANT_SUDO=yes --user root carlosjesusgh/vina_docking:v$(date '+%y%m%d')
docker run --rm -p 8888:8888 --name vina_cont carlosjesusgh/vina_docking:v$(date '+%y%m%d')
docker run --rm --name vina_cont -p 8888:8888 -v ~/vina_docking/example_notebooks:/home/jovyan/example_notebooks carlosjesusgh/vina_docking:v$(date '+%y%m%d')

# push new image to server
docker login -u carlosjesusgh

docker image push carlosjesusgh/vina_docking:v$(date '+%y%m%d') && \
docker image tag carlosjesusgh/vina_docking:v$(date '+%y%m%d') carlosjesusgh/vina_docking:latest && \
docker image push carlosjesusgh/vina_docking:latest

# pull image from server
docker pull docker carlosjesusgh/vina_docking:latest