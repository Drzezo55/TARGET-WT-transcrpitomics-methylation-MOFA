FROM rocker/tidyverse:4.4.3

# Install SSH server and other necessary packages
# Install ALL system dependencies including those needed for Rhtslib
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libcurl4-openssl-dev \
        libz-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        libssl-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libpcre2-dev \
        r-base-dev \
        gcc \
        g++ \
        make \
        cmake \
        git \
        curl \
        zlib1g-dev \
        build-essential \
        openssh-server \
        sudo \
        curl \
        git \
        vim \
        nano \
        htop \
        && \
    # Install renv and devtools
    R -e "install.packages(c('renv', 'devtools'), dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    # Clean up
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/*
 

# Set up SSH server configuration
RUN mkdir -p /var/run/sshd && \
    echo "PasswordAuthentication no" >> /etc/ssh/sshd_config && \
    echo "PermitRootLogin no" >> /etc/ssh/sshd_config && \
    echo "PubkeyAuthentication yes" >> /etc/ssh/sshd_config && \
    echo "X11Forwarding yes" >> /etc/ssh/sshd_config

# Configure sudo for rstudio user (password-less)
RUN echo "rstudio ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# Set up SSH directory for rstudio user
RUN mkdir -p /home/rstudio/.ssh && \
    chown rstudio:rstudio /home/rstudio/.ssh && \
    chmod 700 /home/rstudio/.ssh

# Create renv cache directory
RUN mkdir -p /renv/cache && \
    chown rstudio:rstudio /renv/cache

# Create Positron config directory
RUN mkdir -p /home/rstudio/.positron && \
    chown rstudio:rstudio /home/rstudio/.positron

# Create VS Code config directory for extension compatibility
RUN mkdir -p /home/rstudio/.vscode && \
    chown rstudio:rstudio /home/rstudio/.vscode

# Set up renv configuration
RUN echo 'RENV_PATHS_CACHE=/renv/cache' >> /home/rstudio/.Renviron && \
    echo 'RENV_CONFIG_AUTO_SNAPSHOT=FALSE' >> /home/rstudio/.Renviron && \
    chown rstudio:rstudio /home/rstudio/.Renviron

# Install latest renv
RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

# Expose ports
EXPOSE 22 8787

# Start SSH service and RStudio Server
CMD service ssh start && /init