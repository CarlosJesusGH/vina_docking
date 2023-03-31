#! /bin/bash

alias \
    dockershow="echo '----------------------------------------------------------------------------------------------------' && \
    docker image ls -a && \
    echo '----------------------------------------------------------------------------------------------------' && \
    docker container ls -a && \
    echo '----------------------------------------------------------------------------------------------------' && \
    docker system df && \
    echo '----------------------------------------------------------------------------------------------------' && \
    df -hT && \
    echo '----------------------------------------------------------------------------------------------------' "

alias dockerbuild="docker build . -f Dockerfile -t vina_img:v202303 --progress=plain"
# 

alias \
    dockerrmall="docker container prune -f && \
    docker image prune -f && \
    docker system prune -f && \
    docker builder prune -f && \
    docker system df"

alias dockerrun="docker run -itd --name vina_cont vina_img:v202303"
# alias dockerrun_s02="docker run -itd --name gc3_cont_s02 -p 8000:8000 gc3_img:step02"

# alias dockerterminal="docker exec -it cont_phd_tools bash"
dockerterminal() {
    #do things with parameters like $1 such as
    # mv "$1" "$1.bak"
    # cp "$2" "$1"
    docker exec -it $1 bash
}