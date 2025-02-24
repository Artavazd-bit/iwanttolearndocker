FROM rocker/rbase:latest
FROM rocker/rstudio:latest

RUN apt-get update && apt-get install -y \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libgit2-dev \
    git  \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('renv', 'devtools', 'tidyverse', 'knitr', 'rmarkdown', 'lavaan', 'semTools', 'cSEM', 'stringr'))"

RUN R -e "renv::restore()"





