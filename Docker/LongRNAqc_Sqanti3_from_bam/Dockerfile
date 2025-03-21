# Use the base image with mamba and Sqanti3 already set up
FROM us-central1-docker.pkg.dev/methods-dev-lab/lrtools-sqanti3/lrtools-sqanti3

# Install samtools and BBmap
# RUN mamba install -y samtools

WORKDIR $SRC
RUN wget https://phoenixnap.dl.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && \
     tar -xvzf BBMap_39.01.tar.gz
ENV PATH=$SRC/bbmap/:${PATH}

# Clone the CTAT-LR-fusion repo and the "util" dir to the path
RUN git clone https://github.com/TrinityCTAT/CTAT-LR-fusion.git /CTAT-LR-fusion
ENV PATH="/CTAT-LR-fusion/util:$PATH"


# Create a directory for scripts and add it to PATH
RUN mkdir /scripts
ENV PATH="/scripts:$PATH"

# Copy some custom scripts into the image
COPY scripts/split_bam_per_chromosome.sh /scripts/split_bam_per_chromosome.sh
COPY scripts/split_gtf_per_chromosome.sh /scripts/split_gtf_per_chromosome.sh
COPY scripts/convert_SAM_to_GTF_for_SQANTI3.py /scripts/convert_SAM_to_GTF_for_SQANTI3.py
COPY scripts/concate_gtfs_and_tag_duplicates.py /scripts/concate_gtfs_and_tag_duplicates.py

# Make the scripts excutable
RUN chmod +x /scripts/split_bam_per_chromosome.sh /scripts/split_gtf_per_chromosome.sh /scripts/convert_SAM_to_GTF_for_SQANTI3.py /scripts/concate_gtfs_and_tag_duplicates.py

