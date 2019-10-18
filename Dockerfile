# Base Image
FROM continuumio/miniconda3:4.7.10
 
################## METADATA ######################
 
LABEL base.image="miniconda3:4.7.10"
LABEL version="1"
LABEL software="diag-pipelines-python-r"
LABEL software.version="1.0"
LABEL description="combined python and r packages for diag pipelines"
LABEL tags="Genomics"
 
################## MAINTAINER ######################
 
MAINTAINER Trestan Pillonel
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

RUN conda config --add channels bioconda
RUN conda config --add channels plotly

RUN conda update conda

COPY python-r.yml ./
RUN conda env create -f env.yml
RUN conda clean --all --yes

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /data/
ENV PATH /opt/conda/envs/mummer2circos/bin:$PATH
ENV PATH /usr/local/bin/mummer2circos/:$PATH
