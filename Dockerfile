FROM r-base:4.1.0

# Updates of the system and installation of additional prerequisites.
RUN apt-get update && apt-get install -y \
  file \
  python3-pip \
  virtualenv \
  curl \
  rsync \
  libcurl4-openssl-dev \
  libssl-dev \
  libgtextutils-dev \
# IgBLAST dependencies
  libxml2 \
  libxml2-dev \
  libuv1

# Should we install NCBI-BLAST?
ENV BLAST_INSTALL='Y'

# Main environmental variables
ENV DEPLOYMENT=/programs/deployment
ENV DEPS=$DEPLOYMENT/deps
ENV ARCH=$DEPS/arch
ENV DATA=$DEPLOYMENT/data
ENV IMPORTS=$DEPLOYMENT/imports
ENV PATH=$PATH:$DEPS/flash:$DEPS/cutadapt/bin:$DEPS/fastx/bin:$DEPS/igblast/bin
ENV PATH=$PATH:$DEPS/blast/bin
ENV TEMP=/tmp/deployment

# Stash downloadable archives in a way compatible with docker caching
RUN mkdir -p $TEMP/flash $TEMP/igblast
WORKDIR /archives
RUN curl -SL https://ftp.ncbi.nih.gov/blast/executables/igblast/release/1.18.0/ncbi-igblast-1.18.0-x64-linux.tar.gz \
  | tar -xzC $TEMP/igblast --strip-components 1
RUN curl -SL http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11.tar.gz \
  | tar -xzC $TEMP/flash --strip-components 1

# Optional NCBI-Blast archive
RUN if [ $BLAST_INSTALL = 'Y' ]; then \
  mkdir -p $TEMP/blast; \
  curl -SL https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.11.0/ncbi-blast-2.11.0+-x64-linux.tar.gz \
  | tar -xzC $TEMP/blast --strip-components 1; fi

COPY . /programs

# Set up the deployment environment
RUN mkdir -p $ARCH $DATA $IMPORTS \
             $DEPS/cutadapt $DEPS/fastx/src

# Home system doesn't allow github
# RUN https://github.com/agordon/fastx_toolkit/releases/download/0.0.14/fastx_toolkit-0.0.14.tar.bz2 \
#   | tar -xjC $TEMP/fastx/src
RUN tar xjf $ARCH/fastx_toolkit-0.0.14.tar.bz2 -C $DEPS/fastx/src

# move the stashed executables into final locations
RUN mv $TEMP/*blast $TEMP/flash $DEPS
# RUN mv $TEMP/fastx/* $DEPS/fastx/src

# execute the remaining setup
WORKDIR $DEPLOYMENT
RUN bash ngs-ig_deploy.sh

# clean-up
RUN echo "### Removing the archives." && rm -rv $ARCH/* && \
    echo "### Checking the locations of new binaries ..." && \
    which flash cutadapt fastx_collapser igblastn blastn R && \
    echo "### Done."

CMD bash test.sh
