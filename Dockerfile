FROM debian:stretch

RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bison \
  build-essential \
  cmake \
  flex \
  libgsl0-dev \
  libncurses-dev \
  python-dev \
  python-pip \
  python-setuptools \
  scons \
  libtool \
  libyaml-cpp-dev \
  libz-dev
RUN pip install colored-traceback dendropy matplotlib nestly numpy psutil pysam scipy

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
