# Filename: Dockerfile
# Author: brian.sharber@vumc.org
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/mitochondrial:4.1.0 - deprecated.
# - briansha/maester:4.1.0

FROM ubuntu:18.04

WORKDIR /cromwell_root/

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl \
      zip \
      git \
      unzip \
      gzip \
      g++ \
      make \
      gfortran \
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
      libxml2-dev \
      libssl-dev \
      libfontconfig1-dev \
      libgdal-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libcairo2-dev \
      libxt-dev \
      gdebi-core \
      ca-certificates

# R 4.1.0
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-4.1.0_1_amd64.deb
RUN DEBIAN_FRONTEND=noninteractive gdebi -n r-4.1.0_1_amd64.deb
RUN ln -s /opt/R/4.1.0/bin/R /usr/local/bin/R
RUN ln -s /opt/R/4.1.0/bin/Rscript /usr/local/bin/Rscript

# HTSLib
ADD https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 .
RUN tar -xvf htslib-1.12.tar.bz2 && rm htslib-1.12.tar.bz2
WORKDIR /cromwell_root/htslib-1.12
RUN ./configure --prefix=/usr/local
RUN make
RUN make install
WORKDIR /cromwell_root/
RUN rm -r htslib-1.12

# R packages using a script
#COPY mitochondrial.R .
#RUN Rscript mitochondrial.R

# R packages
RUN R -e "install.packages('tidyverse', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('Matrix', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ggforce', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('data.table', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('readxl', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('Seurat', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('ggrastr', dependencies=TRUE, repos = 'http://cran.us.r-project.org')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'http://cran.us.r-project.org'); library(BiocManager); BiocManager::install('ShortRead')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'http://cran.us.r-project.org'); library(BiocManager); BiocManager::install('SummarizedExperiment')"
RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'http://cran.us.r-project.org'); library(BiocManager); BiocManager::install('ComplexHeatmap')"
