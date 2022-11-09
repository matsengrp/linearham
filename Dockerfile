FROM quay.io/jitesoft/debian:stretch

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
  curl \
  git
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.450-1_amd64.deb && dpkg -i mafft_7.450-1_amd64.deb
RUN Rscript --slave --vanilla -e 'install.packages(c("phylotate", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")'

COPY . /linearham
WORKDIR /linearham

RUN pip install wheel
RUN pip install -r requirements.txt
RUN Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")'

# build from scratch
# RUN cd lib/revbayes/projects/cmake && ./build.sh

# get pre-compiled exec, currently this has compiler issues with the base debian image.
# but the image builds on quay's servers. 

# current error: when you build locally, and 
#lib/revbayes/projects/cmake/rb: /lib/x86_64-linux-gnu/libm.so.6: version `GLIBC_2.29' not found (required by lib/revbayes/projects/cmake/rb)
#lib/revbayes/projects/cmake/rb: /usr/lib/x86_64-linux-gnu/libstdc++.so.6: version `GLIBCXX_3.4.26' not found (required by lib/revbayes/projects/cmake/rb)
#scons: *** [output/cluster-0/lineage_KC576081.1/mcmciter25_mcmcthin1_tuneiter0_tunethin100_numrates4_rngseed0/revbayes_run.trees] Error 1
#scons: building terminated because of errors.

# this happens during the command
# 
RUN curl -fksSL https://github.com/revbayes/revbayes/releases/download/v1.2.1/revbayes-v1.2.1-linux64.tar.gz \
    --output revbayes-v1.2.1-linux64.tar.gz \
    && tar -xvf revbayes-v1.2.1-linux64.tar.gz \
    && (cd lib/revbayes/projects/cmake/ && ln -s /linearham/revbayes-v1.2.1/bin/rb ./)

RUN scons --build-partis-linearham && ./clean.sh

CMD ./test.sh && ./clean.sh
