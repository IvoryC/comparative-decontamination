# hash:sha256:8f3aeefd564dafc045e7acb0c2f3d615aaaeaa5673b90e33714a12e89a78d8c8
FROM registry.codeocean.com/codeocean/r-studio:2022.07.0-548-r4.2.1-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN Rscript -e 'remotes::install_version("glmnet", "4.1.7")' \
    && Rscript -e 'remotes::install_version("splitstackshape", "1.4.8")' \
    && Rscript -e 'remotes::install_version("tidyverse", "2.0.0")' \
    && Rscript -e 'remotes::install_version("torch", "0.10.0")' \
    && Rscript -e 'remotes::install_version("vegan", "2.6-4")'

RUN Rscript -e 'remotes::install_github( \
        "shenhav-and-korem-labs/SCRuB")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "decontam" \
    ))' # Original versions: latest
