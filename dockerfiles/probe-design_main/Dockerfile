# Use a miniforge base
FROM condaforge/miniforge3:23.3.1-0

#update conda channel with bioconda and conda-forge. install software and clean up to reduce size of docker image
RUN conda config --add channels conda-forge && conda config --add channels bioconda && \
conda install -c bioconda bedtools=2.31.1 && \
conda install -c bioconda picard=3.4.0 && \ 
conda install -c bioconda seqkit=2.12.0 && \
conda install -c bioconda bbmap=38.93 && conda clean -afy
