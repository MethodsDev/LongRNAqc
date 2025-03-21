FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $BIN

RUN mamba create --name QC python=3.8
RUN mamba install -n QC numpy pandas seaborn matplotlib -c conda-forge
RUN mamba install -n QC -c bioconda pysam


# Define environment variables for Miniconda to prevent it from asking for confirmation
ENV CONDA_PREFIX=/opt/conda/envs/QC
ENV PATH=$CONDA_PREFIX/bin:$BIN:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false
ENV CONDA_DEFAULT_ENV=QC

# Make RUN commands use the new environment:
SHELL ["mamba", "run", "-n", "QC", "/bin/bash", "-c"]

RUN mamba init

RUN echo "source activate QC" > ~/.bashrc

# Copy script into /usr/local/bin
COPY scripts/violin_plot_indel_rates.py .

# Make script executable
RUN chmod +x violin_plot_indel_rates.py
