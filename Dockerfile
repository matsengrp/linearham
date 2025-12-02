# Multi-stage Dockerfile for linearham
# Stage 1: Extract pre-built RevBayes (with bundled Boost) from base image
# Stage 2: Use modern Debian for Python packages with pre-built wheels
#
# Platform is set to linux/amd64 because partis-bcr contains x86-specific code (SSE2 intrinsics)
# RevBayes with bundled Boost gives ~180x faster MCMC performance than v1.3.0+ without Boost

# Stage 1: Get pre-compiled RevBayes with bundled Boost from Buster base image
FROM --platform=linux/amd64 quay.io/matsengrp/linearham:2025-12-01-base-image AS revbayes-builder

# Stage 2: Modern Debian Bookworm for Python packages
FROM --platform=linux/amd64 debian:bookworm-slim

# Metadata labels
LABEL org.opencontainers.image.source="https://github.com/matsengrp/linearham"
LABEL org.opencontainers.image.description="Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis"
LABEL org.opencontainers.image.licenses="GPL-3.0"
LABEL org.opencontainers.image.title="linearham"

# Set TERM to enable colored-traceback to work in non-TTY environments (e.g., CI)
ENV TERM=xterm

# Copy pre-built RevBayes binary with bundled Boost from the base image
COPY --from=revbayes-builder /usr/local/bin/rb /usr/local/bin/rb

# Install system dependencies needed for linearham
RUN apt-get update && apt-get install -y --no-install-recommends \
  autoconf \
  automake \
  bison \
  build-essential \
  cmake \
  flex \
  gfortran \
  libblas-dev \
  libbz2-dev \
  libgsl-dev \
  liblapack-dev \
  liblzma-dev \
  libtool \
  libyaml-cpp-dev \
  libyaml-dev \
  libz-dev \
  mafft \
  python3-dev \
  python3-pip \
  python3-setuptools \
  r-cran-ape \
  r-cran-coda \
  r-cran-data.table \
  scons \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install R packages needed for phylomd
RUN Rscript --slave --vanilla -e 'install.packages(c("phylotate", "Rcpp", "RcppArmadillo"), repos = "https://cloud.r-project.org")'

# Copy linearham source code
COPY . /linearham
WORKDIR /linearham

# Install Python packages and build linearham
# Use --break-system-packages for Debian Bookworm (externally-managed-environment)
RUN pip3 install --break-system-packages wheel && \
  pip3 install --break-system-packages -r requirements.txt && \
  Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")' && \
  scons --build && ./clean.sh

# Verify RevBayes installation (check if rb is in PATH from the base image)
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD rb --version || exit 1

CMD ./test.sh && ./clean.sh
