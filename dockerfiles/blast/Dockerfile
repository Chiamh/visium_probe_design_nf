# Use a miniforge base
FROM condaforge/miniforge3:23.3.1-0

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda blast=2.12.0 && conda clean -afy
