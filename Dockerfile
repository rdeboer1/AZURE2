FROM --platform=linux/amd64 docker.io/continuumio/miniconda3:latest

LABEL maintainer.name="Jakub Skowronski"
LABEL maintainer.email="jakub.skowronski@pd.infn.it"

WORKDIR /azure2

RUN conda install azure2 -c conda-forge -c skowrons94