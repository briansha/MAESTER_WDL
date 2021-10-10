# Filename: Dockerfile
# Author: brian.sharber@vumc.org
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/maester_homer:4.11

FROM ubuntu:18.04

WORKDIR /HOMER/

RUN apt-get update && apt-get install -y --no-install-recommends \
      curl \
      wget \
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
      gdebi-core \
      ca-certificates

# HOMER - http://homer.ucsd.edu/homer/download.html
RUN curl -O http://homer.ucsd.edu/homer/configureHomer.pl
RUN perl configureHomer.pl -install
RUN cp -r bin /usr/local
