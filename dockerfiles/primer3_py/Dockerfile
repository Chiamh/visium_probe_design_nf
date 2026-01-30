# Use a miniforge base
FROM condaforge/miniforge3:23.3.1-0

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install bioconda::primer3-py=2.3.0 && conda install conda-forge::biopython=1.86 && conda install bioconda::bbmap=39.52 && conda clean -afy
