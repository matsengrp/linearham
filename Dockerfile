FROM quay.io/matsengrp/linearham:2022-11-29-base-image

RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bison \
  libblas-dev \
  build-essential \
  cmake \
  flex \
  gfortran \
  ghostscript \
  graphviz \
  libgsl0-dev \
  liblapack-dev \
  libncurses-dev \
  python-dev \
  python-pip \
  python-setuptools \
  r-cran-ape \
  r-cran-coda \
  r-cran-data.table \
  r-cran-littler \
  scons \
  libtool \
  libyaml-dev \
  libyaml-cpp-dev \
  libz-dev \
  libbz2-dev \
  liblzma-dev \
  less \
  wget \
  git

RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.450-1_amd64.deb && dpkg -i mafft_7.450-1_amd64.deb
RUN Rscript --slave --vanilla -e 'install.packages(c("phylotate", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")'

COPY . /linearham
WORKDIR /linearham

RUN pip install wheel
RUN pip install -r requirements.txt
RUN Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")'
# RUN cd lib/revbayes/projects/cmake && ./build.sh
RUN scons --build-partis-linearham && ./clean.sh

CMD ./test.sh && ./clean.sh
