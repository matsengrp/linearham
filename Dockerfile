# Platform must be explicitly set to linux/amd64 because partis-bcr contains
# x86-specific code (SSE2 intrinsics) that won't compile on ARM64 architectures.
# This ensures consistent builds across different host systems (e.g., Apple Silicon).
ARG TARGETPLATFORM=linux/amd64
FROM --platform=${TARGETPLATFORM} debian:latest

# Set TERM to enable colored-traceback to work in non-TTY environments (e.g., CI)
ENV TERM=xterm

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
  ncurses-base \
  ncurses-term \
  python3-dev \
  python3-pip \
  python3-setuptools \
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
  ca-certificates \
  mafft \
  git
RUN Rscript --slave --vanilla -e 'install.packages(c("phylotate", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")'

# Install RevBayes
RUN mkdir -p /linearham/lib/revbayes/projects/cmake
RUN curl -fksSL https://github.com/revbayes/revbayes/releases/download/v1.2.1/revbayes-v1.2.1-linux64.tar.gz \
    --output /tmp/revbayes-v1.2.1-linux64.tar.gz \
    && tar -xvf /tmp/revbayes-v1.2.1-linux64.tar.gz -C /linearham \
    && rm /tmp/revbayes-v1.2.1-linux64.tar.gz \
    && cd /linearham/lib/revbayes/projects/cmake/ \
    && ln -s /linearham/revbayes-v1.2.1/bin/rb ./

COPY . /linearham
WORKDIR /linearham

RUN pip3 install --break-system-packages wheel
RUN pip3 install --break-system-packages -r requirements.txt
RUN Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")'
RUN scons --build && ./clean.sh

CMD ./test.sh && ./clean.sh
