#!/bin/bash

ENV_PATH="/scratch/drn2/software/conda-fetcher"

# Create environment
conda create -y -p $ENV_PATH python=3.12

# Activate environment
source activate $ENV_PATH

# Core dependencies
conda install -y -c conda-forge \
    requests \
    pandas \
    tqdm \
    numpy \
    scipy \
    xmltramp2 \
    pyyaml \
    networkx

# Bioinformatics packages
conda install -y -c bioconda \
    biopython \
    entrez-direct \
    bioservices \
    biotite

# Database connectors
conda install -y -c conda-forge \
    psycopg2 \
    sqlalchemy \
    pymongo \
    redis-py

# Data processing
conda install -y -c conda-forge \
    h5py \
    pyarrow \
    fastparquet

# Write environment to file
conda env export -p $ENV_PATH > fetcher_environment.yml
