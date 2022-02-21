#!/bin/bash
#
# Initiate an automated run of the pipeline
#
##############################################

DIR="$( cd "$( dirname "$0" )" && pwd )"

nohup bash $DIR/ngs-ig_process.sh process >> $DIR/../run.log
