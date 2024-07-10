#!/bin/bash
#
# Initiate an automated run of the pipeline
#
##############################################

DIR="$( cd "$( dirname "$0" )" && pwd )"

bash $DIR/ngs-ig_process.sh process >> \
    >(tee -a $DIR/../run.log) 2> >(tee -a $DIR/../error.log >&2)
