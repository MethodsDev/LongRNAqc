FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $SRC

# SQANTI3
RUN wget https://github.com/ConesaLab/SQANTI3/archive/refs/tags/v5.3.6.tar.gz && \
     tar -xvf v5.3.6.tar.gz && \
     cd SQANTI3-5.3.6 && \
     sed -i 's/bcbiogff/bcbio-gff/g' SQANTI3.conda_env.yml && \
     mamba env create -f SQANTI3.conda_env.yml -n SQ3



ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -qq update && apt-get -qq -y install \
      automake \
    build-essential \
    bzip2 \
    cmake \
    curl \
    default-jre \
    fort77 \
    ftp \
    g++ \
    gcc \
    gfortran \
    git \
    libblas-dev \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libdb-dev \
    libghc-zlib-dev \
    libjpeg-dev \
    liblzma-dev \
    libncurses-dev \
    libncurses5-dev \
    libpcre3-dev \
    libpng-dev \
    libreadline-dev \
    libreadline-dev \
    libssl-dev \
    libtbb-dev \
    libx11-dev \
    libxml2-dev \
    libxt-dev \
    libzmq3-dev \
    make \
    nano \
    perl \
    pkg-config \
    python3 \
    python3-dev \
    python3-distutils \
    python3-pip \
    python3-setuptools \
    rsync \
    texlive-latex-base \
    tzdata \
    unzip \
    wget \
    x11-common \
    zlib1g-dev





# cDNA Cupcake

SHELL ["mamba", "run", "-n", "SQ3", "/bin/bash", "-c"]

WORKDIR $SRC

RUN git clone https://github.com/MethodsDev/cDNA_Cupcake.git && \
    cd cDNA_Cupcake && \
    python setup.py build && \
    python setup.py install

RUN mamba init

RUN echo "mamba activate SQ3" >> ~/.bashrc

ENV PYTHONPATH=/usr/local/src/cDNA_Cupcake/sequence

ENV PATH=${PATH}:/usr/local/src/SQANTI3-5.3.6

ENV CONDA_DEFAULT_ENV=SQ3
ENV CONDA_PREFIX=/opt/conda/envs/SQ3
ENV PATH=$CONDA_PREFIX/bin:$PATH
