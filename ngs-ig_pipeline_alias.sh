#!/bin/bash
#
#   ngs-ig_pipeline_alias.sh
#   Task: set up the environment for the pipeline script
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# shellcheck disable=SC2034 # Unused variable are invoked in other scripts

RESOURCEDIR="$DEPLOYMENT/imports"

# set the BLAST_INSTALL flag if it has not already been set (possibly, in the Docker container environment)
if [ -z $BLAST_INSTALL ]; then
  BLAST_INSTALL='Y'
fi

## variables passed to FLASH (step 1)
FLASH_minoverlap=20
FLASH_maxoverlap=200
FLASH_mismatch_density=0.45 # 0.25 is the default

## If reads need to be "stitched", this palindromic string is to be used (and recognized later, if needed)
READ_stitch='NNNNNNANNNCNNNGNNNTNNNNNN'

## UMI barcode pattern: "TNNNNTNNNNTNNNNT"
UMIbarcode='T[ATCG]{4}T[ATCG]{4}T[ATCG]{4}T'

## variables passed to cutadapt (step 4)
MINLENGTH=200
MAXLENGTH=800

## IgBLAST variables (step 5)
# It is important for IgBLAST that the IGDATA variable be available in global scope!
export IGDATA="$RESOURCEDIR/igblast_data"
IGBLAST_numthreads=`nproc`

## BLAST variables (optional step used for hinge data)
# BLAST_INSTALL='Y' # expected: Y or N, inheriting variable from Docker container

if [[ $BLAST_INSTALL == 'Y' ]]; then
  export BLAST_DATA="$RESOURCEDIR/blast_data"
  BLAST_numthreads=`nproc`
fi

####################################################################
#  Edit the environmental variables above this line ...
####################################################################
shopt -s expand_aliases

# Check Bash version
if [[ ${BASH_VERSION%%.*} -lt 4 ]]; then
  echo "Error! Bash version >=4 is required. Exiting..."
  exit 1
fi

# Figure out the host system
SYSTEM=`uname`
if [[ $SYSTEM == 'Linux' ]]; then
  grep='grep -P' # enable Perl-like regex
  zcat='zcat'
elif [[ $SYSTEM == 'Darwin' ]]; then
  grep='grep -E'
  zcat='gunzip -c'
else
  echo "Warning! This is an unsupported system. Dragons ahead!"
fi

if [[ ! -e $IGDATA ]] || [[ $BLAST_INSTALL == 'Y' ]]; then
  if [[ ! -e $BLAST_DATA ]]; then
    echo "Error! Check that IgBLAST and BLAST data directories are available by editing the ngs-ig_pipeline_alias.sh ."
    echo "Exiting the pipeline ..."
    exit 1
  fi
fi

## Checking tools
testFlag=0
echo -n "OS................... "; uname -a | cut -d " " -f1
echo -n "Number of CPUs....... "; nproc
echo    "bash................. $BASH_VERSION"
echo -n "perl................. "; perl -v | grep "version"; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "flash................ "; flash -v|head -1; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "fastx_clipper........ "; fastx_clipper -h| grep FASTX; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "fastx_collapser...... "; fastx_collapser -h| grep FASTX; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "fastx_trimmer........ "; fastx_trimmer -h| grep FASTX; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "fastx_quality_stats.. "; fastx_quality_stats -h| grep FASTX; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "fastq_to_fasta....... "; fastq_to_fasta -h| grep FASTX; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "igblastn............. "; igblastn -version| grep Package; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "Rscript.............. "; R --version| head -1; if [[ $? -ne 0 ]]; then testFlag=1; fi
echo -n "cutadapt............. "; cutadapt --version; if [[ $? -ne 0 ]]; then testFlag=1; fi

if [[ $BLAST_INSTALL == 'Y' ]]; then
  echo -n "blastn............... "; blastn -version| grep Package; if [[ $? -ne 0 ]]; then testFlag=1; fi
fi

if [[ $testFlag -ne 0 ]]; then
  echo "Error! Some of the required tools are missing. Please fix the dependencies before continuing."
  echo "Exiting the pipeline ..."
  exit 1
fi
