# Base Image
FROM continuumio/miniconda3:4.10.3
 
################## METADATA ######################
 
LABEL base.image="miniconda3:4.10.3"
LABEL version="1"
LABEL software="mummer2circos"
LABEL software.version="1.1"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Trestan Pillonel
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

RUN conda install -c conda-forge mamba

COPY env.yaml ./
RUN mamba env create -f env.yaml && conda clean --all --yes && echo ok2

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /data/
ENV PATH /opt/conda/envs/mummer2circos/bin:$PATH