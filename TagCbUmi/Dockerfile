# Filename: Dockerfile
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/maester_samtools:1.13

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

# samtools
ADD https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 .
RUN tar -xvf samtools-1.13.tar.bz2 && rm samtools-1.13.tar.bz2
WORKDIR /cromwell_root/samtools-1.13
RUN ./configure --prefix=/usr/local
RUN make
RUN make install
WORKDIR /cromwell_root/
RUN rm -r samtools-1.13
