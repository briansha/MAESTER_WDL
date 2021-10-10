# Filename: Dockerfile
# Used for https://github.com/vangalenlab/MAESTER-2021
# - briansha/star:2.7.9

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
      ca-certificates

# STAR
RUN git clone https://github.com/alexdobin/STAR.git
WORKDIR /cromwell_root/STAR/source
RUN make STAR
RUN cp STAR /usr/local/bin
WORKDIR /cromwell_root/
