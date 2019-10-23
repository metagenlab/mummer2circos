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

RUN conda install conda=4.7.12

COPY env.yaml ./
RUN conda env create -f env.yaml
RUN conda clean --all --yes

RUN mkdir /usr/local/bin/mummer2circos

WORKDIR /usr/local/bin/mummer2circos

COPY mummer2circos.py ./

WORKDIR /usr/local/bin/

RUN git clone https://github.com/metagenlab/TPutils.git

RUN conda init bash
ENTRYPOINT ["/bin/bash"]
WORKDIR /data/
ENV PATH /opt/conda/envs/mummer2circos/bin:$PATH
ENV PATH /usr/local/bin/mummer2circos/:$PATH
ENV PATH /usr/local/bin/TPutils/:$PATH
ENV PYTHONPATH /usr/local/bin/TPutils/:$PYTHONPATH