# Use a full Ubuntu base for legacy Python2 and conda compatibility
FROM ubuntu:20.04

# Use bash for RUN commands
SHELL ["/bin/bash", "-c"]

# Install essential system libraries for old Python / conda envs
RUN apt-get update && apt-get install -y \
    bash \
    bzip2 \
    ca-certificates \
    wget \
    libstdc++6 \
    libgcc-s1 \
    && rm -rf /var/lib/apt/lists/*

# Copy the conda-packed environment from host
COPY probeMining.tar.gz /opt/

WORKDIR /opt

# Extract and fix environment paths
RUN mkdir -p probeMining \
    && tar -xzf probeMining.tar.gz -C probeMining \
    && rm probeMining.tar.gz \
    && /opt/probeMining/bin/python /opt/probeMining/bin/conda-unpack

# Copy OligoMiner scripts
COPY OligoMiner /opt/OligoMiner

# Make scripts executable
RUN chmod +x /opt/OligoMiner/*.py

# Add conda env and OligoMiner scripts to PATH
ENV PATH=/opt/probeMining/bin:/opt/OligoMiner:$PATH

# Optional sanity checks
RUN python --version \
    && echo "Python works"

# Set default workdir (not used by Nextflow for tasks, but useful for interactive runs)
WORKDIR /opt
