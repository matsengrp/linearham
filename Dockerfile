FROM debian:stretch

# RUN apt-get update && apt-get install -y --no-install-recommends \
#   package1 \
#   package2

# RUN conda update -y conda
# RUN conda install -y python=2.7
# RUN conda install -y biopython pandas psutil pysam scons seaborn zlib pyyaml scikit-learn
# RUN conda install -y -c biocore mafft
# RUN pip install colored-traceback dendropy==4.0.0
# COPY . /partis
# WORKDIR /partis
# RUN ./bin/build.sh
# CMD ./test/test.py --quick
