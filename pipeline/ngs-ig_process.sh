#!/bin/bash
#
#   ngs-ig_process.sh
#   Task: main script to run the pipeline
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # # # # # # # # # # # # # # # # # #
#                                   #
#  setting up the environment       #
#                                   #
# # # # # # # # # # # # # # # # # # #

CURDIR="$( cd "$( dirname "$0" )" && pwd )"
WDIR=`dirname $CURDIR`

OUTDIR="00_output"
INDIR="input"
SCRDIR="scripts"
OUT_flash="01_flash_out"
OUT_cutadapt="02_cutadapt_out"
# shellcheck disable=SC2034 # Unused variable are invoked in other scripts
OUT_fastxtk="03_fastx_out"
OUT_igblast="04_igblast_out"

if [[ -f $WDIR/$SCRDIR/ngs-ig_process_functions.sh ]]; then
	# shellcheck source=/dev/null
  source $WDIR/$SCRDIR/ngs-ig_process_functions.sh
else
  echo "Error locating ngs-ig_process_functions.sh"
  exit 1
fi

if [[ -f $WDIR/$SCRDIR/ngs-ig_pipeline_alias.sh ]]; then
  ## import additional environment variables
	# shellcheck source=/dev/null
  source $WDIR/$SCRDIR/ngs-ig_pipeline_alias.sh
else
  echo "Error locating ngs-ig_pipeline_alias.sh"
  exit 1
fi

# User-specific variables collection
if [[ -f $WDIR/$SCRDIR/LabSpecific.sh ]]; then
	# shellcheck source=/dev/null
  source $WDIR/$SCRDIR/LabSpecific.sh
fi

# # # # # # # # # #
#
#  main routine
#
# # # # # # # # # #

## Process commandline and determine whether all processing steps need to be done
case $1 in
  process )
    FRESH=1
    ;;
  repeat )
    echo "###########################"
    echo "# Repeat pipeline run ... #"
    echo "###########################"
    FRESH=0
    ;;
  clean )
    cleanWorkingDirectory
    exit
    ;;
  archive )
    prepareArchive
    exit
    ;;
  *)
    echo "commandline: ngs-ig_process.sh process|repeat|clean|archive"
    exit 0
    ;;
esac

cd $WDIR || { echo "Error: working directory not accessible!"; exit 1; }

## generate SampleManifest.txt with the description of the dataset (from the processing directory name)
if [ ! -f $WDIR/SampleManifest.txt ]; then
  generateProcessManifest
else
  echo "Skipping ... SampleManifest.txt already exists."
fi

## input prep
cd $INDIR || { echo "Error: input directory not accessible!"; exit 1; }
echo "Checking input ..."
if [[ $FRESH -ne 0 ]]; then
  if [[ $? -ne 0 ]]; then
    error "Problem with the input. Check the input directory."
  fi
fi

## capture dataset names
shopt -s nullglob
count=0
for FILE in *.fastq.gz; do
  DATAFILES[$count]=$FILE
  count=$[$count+1]
done

DATA1=${DATAFILES[0]}
DATA2=${DATAFILES[1]}

## determine dataset characteristics from SampleManifest.txt
# shelcheck disable=SC2154
DATANAME=`${grep:?} "Sample manifest for" $WDIR/SampleManifest.txt|sed "s|.* for ||" | sed "s| .*||"`
DATASET_species=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep species | sed "s|.*: ||"`
DATASET_libraryMethod=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep libraryMethod | sed "s|.*: ||"`
DATASET_chain=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep chain | sed "s|.*: ||"`
DATASET_subjectID=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep subjectID | sed "s|.*: ||"`
DATASET_primer=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep primer | sed "s|.*: ||"`
DATASET_libraryType=`$grep -v "^#" $WDIR/SampleManifest.txt | $grep libraryType | sed "s|.*: ||"`

## General warnings about naming conventions
if [[ -z $DATASET_species || -z $DATASET_chain || -z $DATASET_libraryMethod ]]; then
  echo "Error: critical dataset parameters are missing."
  echo "Please make sure that the species/chain/library construction method parameters are specified."
  exit 1
fi

## Sometimes default values are ok.
if [[ -z $DATASET_primer ]]; then
  DATASET_primer="STD"
  echo "#####################################################"
  echo "#                 !!! Warning !!!                   #"
  echo "# The primer entry for this dataset is blank.       #"
  echo "# Please make sure this is not an error.            #"
  echo "# Assuming that this is a 'standard' primer set.    #"
  echo "#####################################################"
fi

if [[ -z $DATASET_libraryType ]]; then
  DATASET_libraryType="variable"
  echo "#####################################################"
  echo "#                  !!! Warning !!!                  #"
  echo "# The library type entry for this dataset is blank. #"
  echo "# Please make sure this is not an error.            #"
  echo "# Assuming that this is a 'variable'-type library.  #"
  echo "#####################################################"
fi

echo "Working on: $DATANAME ..."
echo "subjectID: $DATASET_subjectID"
echo "species: $DATASET_species"
echo "chain: $DATASET_chain"
echo "primer: $DATASET_primer"
echo "libraryMethod: $DATASET_libraryMethod"
echo "libraryType: $DATASET_libraryType"

###############################################
## initialize the output directory structure ##
###############################################
if [[ $FRESH -ne 0 ]]; then
  if [[ -f $WDIR/running.txt ]]; then
    echo "Already running a process! This was probably started by mistake."
    exit 1
  fi
  initialize
fi

checkExist "$WDIR/$SCRDIR/adapter*.conf"
if [[ $? -ne 1 ]]; then
  selectAdaptors $DATASET_species $DATASET_chain $DATASET_libraryMethod $DATASET_primer $DATASET_libraryType
else
  echo "Skipping ... Adapters selected during previous run."
fi

##################################
## generate FASTQ quality plots ##
##################################
checkExist "$WDIR/$OUTDIR/*.quality.pdf"
if [[ $? -ne 1 ]]; then
  plotRunQuality $DATA1 $DATA2
else
  echo "Skipping FASTQ quality plot generation."
fi

##################
#### FLASH step ##
##################
# checkExist "$WDIR/$OUT_flash/out.extendedFrags.fastq*"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUT_flash/flash_stage_complete.txt ]]; then
  echo "Skipping step ... amplicons reconstructed during previous run."
else
  FLASHstep ${FLASH_maxoverlap:?} ${FLASH_minoverlap:?} ${FLASH_mismatch_density:?} \
    "$WDIR/$INDIR/$DATA1" "$WDIR/$INDIR/$DATA2"
  touch $WDIR/$OUT_flash/flash_stage_complete.txt
fi

###################
## CUTADAPT step ##
###################
# checkTargetNewer "$WDIR/$OUT_flash/out.extendedFrags.fastq*" "$WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq*"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUT_cutadapt/cutadapt_stage_complete.txt ]]; then
  echo "Skipping step ... Primers trimmed during previous run."
else
  cutadaptStep
  touch $WDIR/$OUT_cutadapt/cutadapt_stage_complete.txt
fi

########################
## fastx-toolkit step ##
########################
# checkTargetNewer "$WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq*" "$WDIR/$OUT_igblast/input*.fasta"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUT_fastxtk/fastxtk_stage_complete.txt ]]; then
echo "Skipping step ... Dataset converted and collapsed during previous run."
else
  fastxStep $DATASET_libraryMethod $DATASET_primer $DATASET_libraryType
  touch $WDIR/$OUT_fastxtk/fastxtk_stage_complete.txt  
fi

###################
## igblastn step ##
###################
# checkTargetNewer "$WDIR/$OUT_igblast/input*.fasta" "$WDIR/$OUT_igblast/$DATANAME.aa.*igblast_out*"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUT_igblast/igblast_stage_complete.txt ]]; then
  echo "Skipping step ... IgBLAST annotation completed during previous run."
else
  IgBLASTstep $DATASET_species
  touch $WDIR/$OUT_igblast/igblast_stage_complete.txt
fi

###############################
## IgBLAST output processing ##
###############################
# checkTargetNewer "$WDIR/$OUT_igblast/$DATANAME.aa.*igblast_out*" "$WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta*"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUT_igblast/igblast_output_processing_stage_complete.txt ]]; then
  echo "Skipping step ... IgBLAST output processing completed during previous run."
else
  if [ ! -f $WDIR/$OUT_igblast/${DATANAME}.aa.igblast_out.airr.tsv ] && [ ! -f $WDIR/$OUT_igblast/${DATANAME}.aa.read1.igblast_out.airr.tsv ]; then
    gunzip $WDIR/$OUT_igblast/${DATANAME}.*.igblast_out.airr.tsv.gz
    if [[ $? -ne 0 ]]; then
      error "could not locate the IgBLAST output for repeat processing."
    fi
  fi
  IgBLASToutputProcessing $DATASET_chain $DATASET_primer $DATASET_libraryMethod $DATASET_libraryType
  touch $WDIR/$OUT_igblast/igblast_output_processing_stage_complete.txt
fi

###########################
## Hinge processing step ##
###########################
if [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano)$ ]] ; then
  # checkTargetNewer "$WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta*" \
  #   "$WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.subclass.subset.CDR3aa_dict.fasta"
  # TODO: this doesn't take care of the newer UMI5RACE datasets, only the multiplex ones
  # if [[ $? -ne 1 ]]; then
    if [[ -f $WDIR/$OUT_igblast/hinge_processing_stage_complete.txt ]]; then
    echo "Skipping step ... hinge BLAST output processing completed during previous run."
  else
    hingeProcessingStep $DATASET_species
    touch $WDIR/$OUT_igblast/hinge_processing_stage_complete.txt
  fi
fi

#######################################
## Post-processing and visualization ##
#######################################
# checkTargetNewer "$WDIR/$OUTDIR/$DATANAME.process_stats" "$WDIR/$OUTDIR/FASTAViewer/*.RData"
# if [[ $? -ne 1 ]]; then
if [[ -f $WDIR/$OUTDIR/postprocessing_stage_complete.txt ]]; then
  echo "Skipping post-processing steps (figures and FASTAViewer processing)."
else
  if [ -f $WDIR/$SCRDIR/postprocess/postprocess.sh ]; then
    echo "########################################"
    echo "Post-processing (visualization, etc.) ..."
		# shellcheck source=/dev/null
    source $WDIR/$SCRDIR/postprocess/postprocess.sh
    touch $WDIR/$OUTDIR/postprocessing_stage_complete.txt
  fi
fi

############################
## Compress intermediates ##
############################
compressIntermediates

if [[ $FRESH -ne 0 ]]; then
  mv $WDIR/running.txt $WDIR/done.txt
  echo "done" >> $WDIR/done.txt
  date >> $WDIR/done.txt
  cat $WDIR/done.txt >> $WDIR/run.log
fi
