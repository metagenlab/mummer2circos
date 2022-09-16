# Base Image
FROM mambaorg/micromamba:0.24.0
 
################## METADATA ######################
 
LABEL base.image="miniconda3:4.12.0"
LABEL software="mummer2circos"
LABEL software.version="1.4.2"
LABEL tags="Genomics"
 
################## INSTALLATION ######################
ENV DEBIAN_FRONTEND noninteractive

RUN micromamba install --yes --name base --channel conda-forge --channel bioconda mummer2circos && micromamba clean --all --yes

#RUN conda init bash
#ENTRYPOINT ["/bin/bash"]
WORKDIR /data/
ENV PATH /opt/conda/bin:$PATH