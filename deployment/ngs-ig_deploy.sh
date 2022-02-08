#!/bin/bash
# NGS-Ig deployment script
#

# Build section
echo "### Building FLASH..."
cd $DEPS/flash || { echo "Error: Could not find FLASH source directory!"; exit 1; }
make |tee -a $DEPS/FLASH_make.log
echo "### Done."

echo "### Building FASTX-toolkit..."
cd $DEPS/fastx/src/fastx_toolkit-0.0.14 || { echo "Error: Could not find fastx source directory!"; exit 1; }
patch $DEPS/fastx/src/fastx_toolkit-0.0.14/src/fasta_formatter/fasta_formatter.cpp $DEPLOYMENT/fastx-toolkit-gcc7-patch.txt
patch $DEPS/fastx/src/fastx_toolkit-0.0.14/src/libfastx/fastx.h $DEPLOYMENT/fastx_patch.txt
./configure --prefix=$DEPS/fastx |tee -a $DEPS/fastx_configure.log
make |tee -a $DEPS/fastx_make.log
make install |tee -a $DEPS/fastx_install.log
echo "### Done."

# Python setup section
echo "### Installing cutadapt using pip ..."
cd $DEPS/cutadapt || { echo "Error: Could not find cutadapt directory!"; exit 1; }
virtualenv $DEPS/cutadapt
$DEPS/cutadapt/bin/pip install cutadapt |tee -a $DEPS/cutadapt_install.log
echo "### Done."

# R package setup section
echo  "### Installing R packages ..."
# needed packages:
# dada2
# optparse
# here
# ggplot2

R -e 'install.packages(c("here","optparse","ggplot2","BiocManager","png","jpeg","reshape2","latticeExtra","matrixStats","bitops","Rcpp","RcppParallel","RCurl","hwriter"))' |tee -a $DEPS/Rpackages_install.log
R -e 'BiocManager::install("dada2")'  |tee -a $DEPS/Rpackages_install.log

R CMD INSTALL -l /usr/local/lib/R/site-library $ARCH/Rpackages/* |tee -a $DEPS/Rpackages_install.log
