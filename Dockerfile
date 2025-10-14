# Author: Roberto Olvera Hernandez
# Date: 2025-10-14
FROM rocker/verse:4.3.1

WORKDIR /tmp/

# Install system dependencies (e.g., for 'sf', 'curl', 'xml2' packages)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev && \
    rm -rf /var/lib/apt/lists/*

# Copy the list of R packages to install first (for better caching)
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Copy the rest of the pipeline code
COPY utils/ utils/
COPY gwas-make-plots.R gwas-make-plots.R
COPY map-genes.R map-genes.R

# Install renv and restore the project library
RUN Rscript -e 'install.packages(c("renv", "BiocManager"))'
RUN Rscript -e 'renv::restore(lockfile = "renv.lock", repos = NULL)'