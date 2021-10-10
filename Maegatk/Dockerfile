# Filename: Dockerfile
# Author: brian.sharber@vumc.org
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/maester_maegatk:v01

FROM ubuntu:18.04

WORKDIR /cromwell_root/

# Install Python3, Java, and other packages.
RUN apt-get update && apt-get install -y --no-install-recommends \
      wget \
      curl \
      zip \
      git \
      unzip \
      gzip \
      g++ \
      make \
      zlib1g-dev \
      libgfortran4 \
      liblapacke-dev \
      libopenblas-dev \
      libbz2-dev \
      liblzma-dev \
      libcurl4-openssl-dev \
      libpng-dev \
      libjpeg-dev \
      zlib1g-dev \
      gdebi-core \
      ca-certificates \
      python3 \
      python3-pip \
      python3-setuptools \
      python3-dev \
      default-jre \
      libncurses5-dev

# Maegatk seems to attempt to run some things using the python command, when it actually needs to run python3.
# Alias python3 as python.
# https://stackoverflow.com/questions/36388465/how-to-set-bash-aliases-for-docker-containers-in-dockerfile
# Dockerfile - installation of packages uses a non-interactive shell
# WDL command section - using 'alias python=python3' did not work. 
#   - The symbolic link at the bottom of this script below (ln -s) does.
# Many scripts use the python command often in their installation, such as bedtools using the python command in its make script.
# However, in maegatk's case, it needs to run python3 whenever it wants to run python. (As some packages it wants to use are only available for python3.)
# The RUN command below works, but only for installation of packages onto the docker image.
#   - In the WDL command section, it is prone to some errors:
#     - /usr/bin/python: 1: /usr/bin/python: -e: not found
#     - python3
#RUN echo 'alias python=python3' >> ~/.bashrc
RUN echo -e '#!/bin/bash\necho python3' > /usr/bin/python && \
      chmod +x /usr/bin/python

# bwa
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /cromwell_root/bwa
RUN make
RUN cp bwa /usr/local/bin
WORKDIR /cromwell_root/

# bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
RUN tar -zxvf bedtools-2.29.1.tar.gz
WORKDIR /cromwell_root/bedtools2
RUN make
WORKDIR /cromwell_root/bedtools2/bin
RUN cp * /usr/local/bin
WORKDIR /cromwell_root/

# maegatk - https://github.com/caleblareau/maegatk
RUN pip3 install maegatk

# The two below need to be in the WDL command section for maegatk to run properly.
# export LC_ALL=C.UTF-8
# export LANG=C.UTF-8

# Libraries maegatk needs - pysam, R 4.1.0, BiocManager (SummarizedExperiment), data.table, setuptools, ruamel_yaml, samtools
#   - These are found out not in the documention, but rather when trying to run an analysis...and it fails.
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-4.1.0_1_amd64.deb
RUN DEBIAN_FRONTEND=noninteractive gdebi -n r-4.1.0_1_amd64.deb
RUN ln -s /opt/R/4.1.0/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.1.0/bin/Rscript /usr/local/bin/Rscript
RUN pip3 install pysam
RUN R -e "install.packages('dplyr', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'http://cran.us.r-project.org')"
RUN R -e "BiocManager::install('SummarizedExperiment')"
RUN R -e "install.packages('data.table', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN pip3 install setuptools
RUN pip3 install ruamel_yaml
ADD https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 .
RUN tar -xvf samtools-1.13.tar.bz2 && rm samtools-1.13.tar.bz2
WORKDIR /cromwell_root/samtools-1.13
RUN ./configure --prefix=/usr/local
RUN make
RUN make install
WORKDIR /cromwell_root/
RUN rm -r samtools-1.13

# Remove the python alias - make a symbolic link instead - this will stick.
RUN rm /usr/bin/python
RUN ln -s /usr/bin/python3 /usr/bin/python
