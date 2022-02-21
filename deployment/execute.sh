#!/bin/bash

DIR=$1

if [ ! -d /mnt/$DIR ] ; then
    echo "Error! Cannot find $DIR"
    exit 1
fi

cd /mnt/$DIR || { echo "Error: working directory not accessible!"; exit 1; }

# check if there is a custom pipeline provided in the scripts directory
if [ ! -d ./scripts ] ; then
   echo "No scripts directory provided. Using the ngs-ig pipeline included in the container."
   rsync -av /programs/pipeline/* ./scripts --exclude deployment --exclude .git*
else
   echo "Found the custom scripts included with the input data."
fi


bash scripts/ngs-ig_process.sh process 2>&1 |tee ./run.log

echo "Done with the run. Check the output directory."
