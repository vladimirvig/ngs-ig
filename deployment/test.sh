#!/bin/bash
# this is intended to be a hard-coded test requiring read-write access to the container (i.e., not for singularity)

DATASET="Mm-BALB-p6wk-spleen-UMI5RACENEB-variableNano-vv874pool876-IgG"
DATASET_DATE=20200826
TARGET="/programs/deployment/data/${DATASET_DATE}/${DATASET}"

cd $TARGET || { echo "Error: working directory not accessible!"; exit 1; }
rsync -av /programs/* $TARGET/scripts --exclude deployment --exclude .git*
bash $TARGET/scripts/ngs-ig_process.sh process 2>&1 |tee ./run.log

mkdir -p /mnt/${DATASET_DATE}_${DATASET}_out
cp -ir $TARGET/* /mnt/${DATASET_DATE}_${DATASET}_out
echo "Done with the test run. Check the output summary."
