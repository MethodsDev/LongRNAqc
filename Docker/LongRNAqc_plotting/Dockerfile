FROM condaforge/mambaforge

ENV SRC=/usr/local/src
ENV BIN=/usr/local/bin

WORKDIR $BIN

RUN mamba create --name QC python=3.11
RUN mamba install -n QC numpy pandas seaborn matplotlib colorcet pgzip -c conda-forge


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
COPY scripts/LongRNAqc_classification_plots.py .

# Make script executable
RUN chmod +x LongRNAqc_classification_plots.py
