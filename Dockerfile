# Unified Dockerfile for linearham using official RevBayes Docker image
# Platform is set to linux/amd64 because partis-bcr contains x86-specific code (SSE2 intrinsics)
# Note: This will be slow on Apple Silicon due to emulation. For local development on ARM Macs,
# consider running tests in CI instead, or accept that MCMC tests will be very slow.
ARG TARGETPLATFORM=linux/amd64
FROM --platform=${TARGETPLATFORM} sswiston/phylo_docker:full_amd64

# Metadata labels
LABEL org.opencontainers.image.source="https://github.com/matsengrp/linearham"
LABEL org.opencontainers.image.description="Bayesian phylogenetic hidden Markov model for B cell receptor sequence analysis"
LABEL org.opencontainers.image.licenses="GPL-3.0"
LABEL org.opencontainers.image.title="linearham"

# Set TERM to enable colored-traceback to work in non-TTY environments (e.g., CI)
ENV TERM=xterm

# Install additional system dependencies needed for linearham
RUN apt-get update && apt-get install -y --no-install-recommends \
    autoconf \
    automake \
    bison \
    flex \
    libtool \
    libyaml-dev \
    libyaml-cpp-dev \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    scons \
    mafft \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Copy linearham source code
COPY . /linearham
WORKDIR /linearham

# Install Python packages, R packages, and build linearham in a single layer
RUN pip3 install --break-system-packages wheel && \
    pip3 install --break-system-packages -r requirements.txt && \
    Rscript --slave --vanilla -e 'install.packages("lib/phylomd", repos = NULL, type = "source")' && \
    scons --build && ./clean.sh

# Verify RevBayes installation (check if rb is in PATH from the base image)
HEALTHCHECK --interval=30s --timeout=3s --start-period=5s --retries=3 \
  CMD rb --version || exit 1

CMD ./test.sh && ./clean.sh
