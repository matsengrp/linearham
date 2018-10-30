FROM debian:stretch

RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bison \
  libblas-dev \
  build-essential \
  cmake \
  flex \
  libgsl0-dev \
  libncurses-dev \
  python-dev \
  python-pip \
  python-setuptools \
  r-cran-ape \
  r-cran-data.table \
  r-cran-littler \
  scons \
  libtool \
  libyaml-dev \
  libyaml-cpp-dev \
  libz-dev
RUN pip install colored-traceback dendropy jinja2 matplotlib nestly numpy psutil pysam pyyaml scipy

COPY . /linearham
WORKDIR /linearham


# RUN conda update -y conda
# RUN conda install -y python=2.7
# RUN conda install -y biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn
# RUN conda install -y -c biocore mafft
# RUN pip install colored-traceback dendropy==4.0.0
# COPY . /partis
# WORKDIR /partis
# RUN ./bin/build.sh
# CMD ./test/test.py --quick
