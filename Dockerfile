FROM rocker/shiny:4.0.5

# Install system requirements for index.R as needed
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    git-core \
    libssl-dev \
    libcurl4-gnutls-dev \
    curl \
    libsodium-dev \
    libxml2-dev \
    libicu-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*


####### System libraries #######

RUN sudo apt-get update && apt-get install -y \
    --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \ 
    git-core \
    curl \
    libsodium-dev \
    libxml2-dev \
    libicu-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN apt-get update && apt-get install -y --no-install-recommends \
	libmariadb-dev

RUN sudo apt-get update && apt-get install -y --no-install-recommends \
	libfreetype6-dev \
	libpng-dev \
	libtiff5-dev \
	libjpeg-dev

RUN apt-get update && apt-get install -y --no-install-recommends \
  libudunits2-dev \
  libharfbuzz-dev \
  libfribidi-dev 
  
RUN apt-get update && apt-get install -y --no-install-recommends \
	cmake
	

####### R dependencies #######

RUN R -e 'install.packages(c("AER", "concaveman", "nloptr","units", "pbkrtest", "sf", "car", "terra", "lme4"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shiny", "shinydashboard", "ggstance", "bslib", "tidyverse", "BiocManager"),repos="http://cran.rstudio.com/", dependencies=T)'

RUN R -e 'BiocManager::install("GenomicRanges")'

RUN R -e 'install.packages(c("DBI", "shinyjs", "shinycssloaders", "shinyBS", "shinydashboard"), repos="http://cran.rstudio.com/", dependencies=T)'

RUN R -e 'install.packages(c("stringr", "data.table", "ggforce", "gridExtra", "sandwich"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("DT"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shinylogs"), repos="http://cran.rstudio.com/", dependencies = T)'

######## ggtranscript ##########

RUN R -e 'install.packages(c("devtools"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'devtools::install_github("dzhang32/ggtranscript")'

####### other R libraries ######

RUN R -e 'install.packages(c("pbkrtest"), source="https://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("ggpubr"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("ggrepel"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'install.packages(c("shinybusy"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("tidytext"), repos="http://cran.rstudio.com/", dependencies = T)'
RUN R -e 'install.packages(c("here"), repos="http://cran.rstudio.com/", dependencies = T)'

RUN R -e 'remove.packages("ggplot2")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.0.tar.gz", repos=NULL, type="source")'

RUN R -e 'remove.packages("ggrepel")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/ggrepel/ggrepel_0.9.1.tar.gz", repos=NULL, type="source")'

####### COPY Rprofile #########

COPY Rprofile.site /usr/lib/R/etc/
RUN mkdir /root/introverse_v2
COPY . /root/introverse_v2
RUN mkdir /root/introverse_v2/database

####### EXPOSE #########

EXPOSE 3838


CMD ["/usr/bin/shiny-server"]
