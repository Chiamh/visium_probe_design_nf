# Use a miniforge base
FROM condaforge/miniforge3:23.3.1-0

# Use bash for RUN commands
SHELL ["/bin/bash", "-c"]

# Environment name
ARG ENV_NAME=fasta_download
ENV ENV_NAME=${ENV_NAME}

# Configure conda channels (order matters!)
RUN conda config --system --add channels conda-forge \
	&& conda config --system --add channels bioconda \
	&& conda config --system --set report_errors false

#Disable conda notices/spinner globally
ENV CONDA_NO_PLUGINS=true
ENV CONDA_NOTICES=false

RUN conda create -y -n ${ENV_NAME} \
	r-base \
	r-tidyverse \
	r-stringr \
	r-biomartr \
	r-optparse \
	&& conda clean -afy

# Activate environment by default
ENV PATH=/opt/conda/envs/${ENV_NAME}/bin:$PATH
