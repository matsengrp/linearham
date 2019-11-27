FROM debian:stretch

RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bison \
  libblas-dev \
  build-essential \
  cmake \
  flex \
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
  wget
RUN wget https://mafft.cbrc.jp/alignment/software/mafft_7.450-1_amd64.deb && dpkg -i mafft_7.450-1_amd64.deb
RUN pip install biopython colored-traceback dendropy graphviz jinja2 matplotlib nestly numpy psutil pysam pyyaml scipy weblogo
RUN Rscript --slave --vanilla -e 'install.packages(c("phylotate", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")'

COPY . /linearham
WORKDIR /linearham

RUN Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")'
RUN cd lib/revbayes/projects/cmake && ./build.sh
RUN scons --build-partis-linearham

CMD scons --run-partis --fasta-path=data/liao_dataset.fasta --all-clonal-seqs \
    && scons --run-linearham --template-path=templates/revbayes_template.rev --mcmc-iter=25 --mcmc-thin=1 --tune-iter=0 --lineage-unique-id=KC576081.1
