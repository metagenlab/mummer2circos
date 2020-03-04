# Base Image
FROM continuumio/miniconda3:4.7.10
 
################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.10"
LABEL version="1"
LABEL software="mummer2circos"
LABEL software.version="1.1"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Trestan Pillonel
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

RUN conda install conda=4.7.12

COPY env.yaml ./
RUN conda env create -f env.yaml
RUN conda clean --all --yes

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /data/
ENV PATH /opt/conda/envs/mummer2circos/bin:$PATH