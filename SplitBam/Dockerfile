# Filename: Dockerfile
# Author: brian.sharber@vumc.org
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/bamtools:latest

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
      libncurses5-dev \
      gdebi-core \
      ca-certificates

# bamtools
RUN apt-get update && apt-get install -y --no-install-recommends bamtools
