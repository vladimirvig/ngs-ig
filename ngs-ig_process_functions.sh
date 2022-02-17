#!/bin/bash
#
#   ngs-ig_process_functions.sh
#   Library of functions for pipeline operation
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# function: error handling
# arguments: error message
function error () {
  echo "Error detected: $1"
  echo "Exiting the pipeline..."
  if [ -f $WDIR/running ]; then
    mv -f $WDIR/running $WDIR/error
  fi
  exit 1
}

function time_msg () {
  echo "[$(date +%T%Z)]...$1"
}

# function: clean directory
# arguments: none
function cleanWorkingDirectory () {
  cd $WDIR || { echo "Error: working directory not accessible!"; exit 1; }
  echo "Deleting processing subdirectories ..."
  rm -rf 0*
  echo "Deleting pipeline status files ..."
  rm -f "done" "error" "running"
  echo "Deleting SampleManifest.txt..."
  rm -f SampleManifest.txt
  echo "Deleting run.log ..."
  rm -f run.log
  echo "Deleting archiving.log"
  rm -f archiving.log
  echo "Deleting adapter selections..."
  rm -f $SCRDIR/adapter?.*
  echo "Exiting..."
}

# function: setup
# arguments: none
function initialize () {
  cd $WDIR || { echo "Error: working directory not accessible!"; exit 1; }
  echo "run start:" > $WDIR/running
  date >> $WDIR/running

  mkdir $OUTDIR
  if [[ $? -ne 0 ]]; then
    error "Problem setting up workspace. Check that your directory only contains \"input\" and \"scripts\" directories."
  fi
  mkdir ${OUT_flash:?}
  if [[ $? -ne 0 ]]; then
    error "Problem setting up workspace. Check that your directory only contains \"input\" and \"scripts\" directories."
  fi
  mkdir ${OUT_cutadapt:?}
  if [[ $? -ne 0 ]]; then
    error "Problem setting up workspace. Check that your directory only contains \"input\" and \"scripts\" directories."
  fi
  mkdir ${OUT_fastxtk:?}
  if [[ $? -ne 0 ]]; then
    error "Problem setting up workspace. Check that your directory only contains \"input\" and \"scripts\" directories."
  fi
  mkdir ${OUT_igblast:?}
  if [[ $? -ne 0 ]]; then
    error "Problem setting up workspace. Check that your directory only contains \"input\" and \"scripts\" directories."
  fi
}

# function: checkExist
# arguments: path to list of files
# returns: 1 if found
function checkExist (){
  list=$1
  value=0

  for file in $list; do
    if [[ -f $file ]]; then
      return 1
    fi
  done
  return $value
}

# function: checkTargetNewer
# arguments: path to list of input files, path to list of output (target) files
# returns: 1 if true
function checkTargetNewer (){
  inList=$1
  outList=$2
  value=0

  for inFile in $inList; do
    if [[ -f $inFile ]]; then
      for outFile in $outList; do
        if [[ "$outFile" -nt "$inFile" ]]; then
          return 1
        fi
      done
    fi
  done
  return $value
}

# function: compressIntermediates
# arguments: none
function compressIntermediates () {
  echo "####################################################################"
  time_msg "Compressing the intermediate FASTQ and *.igblast_out files ..."

  find $WDIR -name "*.fastq" -exec gzip -f -9 {} \; -exec echo -n "." \;
  find $WDIR -name "*.*blast_out" -exec gzip -f -9 {} \; -exec echo -n "." \;
  echo ""
}

# function: prepareArchive
# arguments: none
function prepareArchive () {
  compressIntermediates
	time_msg "Compressing the fasta files ..."
  find $WDIR -name "*.fasta" -exec gzip -f -v -9 {} \; 2>&1 | tee -a $WDIR/archiving.log
}

# function: generate SampleManifest.txt
# arguments: none
function generateProcessManifest () {
  bash $WDIR/$SCRDIR/build_sample_manifest.sh $WDIR > $WDIR/SampleManifest.txt
  if [[ $? -ne 0 ]]; then
    echo "Sample manifest problem detected ..."
    cat $WDIR/SampleManifest.txt
    error "Couldn't parse dataset name."
  fi
}

# function: selectAdaptors
# arguments: species, chain, libraryMethod, primer, libraryType
function selectAdaptors (){
  species=$1
  chain=$2
  libraryMethod=$3
  primer=$4
  libraryType=$5

  ## set the variable for cutadapt 5' adapter
  ##       "species:chain:libraryMethod"
  case $libraryMethod in
    UMI5RACE)
      adapter5_lookup="all:all:UMI5RACE"
      ;;
    UMI5RACEASYM)
      adapter5_lookup="all:all:UMI5RACEASYM"
      adapter3_lookup="all:all:UMI5RACEASYM"
      cat $WDIR/$SCRDIR/adapters/primers/adapterUMI5RACE.fasta > $WDIR/$SCRDIR/adapter5.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter${primer}.fasta >> $WDIR/$SCRDIR/adapter5.fasta
      ;;
    UMI5RACENEB)
      adapter5_lookup="all:all:UMI5RACENEB"
      adapter3_lookup="all:all:UMI5RACENEB"
      cat $WDIR/$SCRDIR/adapters/primers/adapterUMI5RACE.fasta > $WDIR/$SCRDIR/adapter5.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter${primer}.fasta >> $WDIR/$SCRDIR/adapter5.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapterUMI5RACE_revcomp.fasta > $WDIR/$SCRDIR/adapter3.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter${primer}_revcomp.fasta >> $WDIR/$SCRDIR/adapter3.fasta
      ;;
    multiplexNEB)
      adapter5_lookup="all:all:UMI5RACENEB"
      adapter3_lookup="all:all:UMI5RACENEB"
      cat $WDIR/$SCRDIR/adapters/primers/adapter5multiplex.fasta > $WDIR/$SCRDIR/adapter5.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter${primer}.fasta >> $WDIR/$SCRDIR/adapter5.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter5multiplex_revcomp.fasta > $WDIR/$SCRDIR/adapter3.fasta
      cat $WDIR/$SCRDIR/adapters/primers/adapter${primer}_revcomp.fasta >> $WDIR/$SCRDIR/adapter3.fasta
      ;;
    *)
      error "unrecognized library prep method $libraryMethod in setup! Check that the dataset name follows correct nomenclature."
  esac

  # validity check for primers/adapters
  if [ ${primers5[${adapter5_lookup}]} ]
  then
    adapter5=${primers5[${adapter5_lookup}]}
  else
    error "invalid primer chosen ... $adapter5_lookup"
  fi

  if [ ${primers3[${adapter3_lookup}]} ]
  then
    adapter3=${primers3[${adapter3_lookup}]}
  else
    error "invalid primer chosen ... $adapter3_lookup"
  fi

  ## copying the adapter files for cutadapt
  time_msg "Setup: Using $adapter5 and $adapter3 as adapters."
  cp $WDIR/$SCRDIR/adapters/$adapter5 $WDIR/$SCRDIR/adapter5.conf
  cp $WDIR/$SCRDIR/adapters/$adapter3 $WDIR/$SCRDIR/adapter3.conf
}

# function: plotRunQuality
# arguments: file1 file2
function plotRunQuality () {
  file1=$1
  file2=$2

  ## Generate sequencing run statistics and graphs
  time_msg "Generating stats for $file1"
  Rscript $WDIR/$SCRDIR/readQualityPlot.R --inputFile=$WDIR/$INDIR/$file1 --outputFile=$WDIR/$OUTDIR/$file1.quality.pdf
  time_msg "Generating stats for $file2"
  Rscript $WDIR/$SCRDIR/readQualityPlot.R --inputFile=$WDIR/$INDIR/$file2 --outputFile=$WDIR/$OUTDIR/$file2.quality.pdf

  ## accounting start
  buildReadAccountingSummary initialize
  buildUniqueAccountingSummary initialize
}

# function: the FLASH step
# arguments: maxoverlap, minoverlap, mismatchDensity, file1, file2
function FLASHstep (){
  maxoverlap=$1
  minoverlap=$2
  mismatchDensity=$3
  file1=$4
  file2=$5

  cd $WDIR/$OUT_flash || { error "Error: FLASH output directory not accessible!"; }
  if [[ "${DATASET_libraryType:?}" =~ ^(variableNano|HINGENano)$ ]] || ([[ "${DATASET_libraryMethod:?}" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" == HINGE ]]); then
    flash -M $maxoverlap -m $minoverlap -x $mismatchDensity -z $file1 $file2 --output-prefix=out1
    time_msg "Generating a stitched FASTQ from the paired reads."
    gunzip -c $file1 > temp1.fastq
    gunzip -c $file2 > temp2.fastq
    perl $WDIR/$SCRDIR/fastq_stitch.pl temp1.fastq temp2.fastq ${READ_stitch:?} > out.extendedFrags.fastq
    gzip out.extendedFrags.fastq
  else
    flash -M $maxoverlap -m $minoverlap -x $mismatchDensity -z $file1 $file2
  fi

  ## accounting
  buildReadAccountingSummary extendStepAcct
  buildUniqueAccountingSummary extendStepAcct
}

# function: the FLASH step for asymmetric sequencing experiment
# arguments: minoverlap, mismatchDensity, interleaved FASTQ file
function asymmetricSequencingExtension (){
  maxoverlap=$1
  minoverlap=$2
  mismatchDensity=$3
  file=$4
  workingDir=`pwd`
  cd $WDIR/$OUT_flash || { error "Error: FLASH output directory not accessible!"; }
  flash -M $maxoverlap -m $minoverlap -x $mismatchDensity --interleaved-input -o UMI5RACEASYM -z $file
  cd $workingDir || { error "Error: working directory no longer accessible!"; }
}

# function: the cutadapt step
# arguments:
function cutadaptStep (){
  cd $WDIR/$OUT_cutadapt || { error "Error: Cutadapt output directory not accessible!"; }
  echo "####################################"
  time_msg "Running cutadapt for 5' end..."
  args=$(tr "\n" " " <$WDIR/$SCRDIR/adapter5.conf)
  if [[ "$DATASET_libraryMethod" == multiplexNEB ]]; then
    if [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
      args="$args --no-trim"
    fi
    echo "Hinge dataset from NEB-adaptored library: preserving the forward (multiplex) primer sequences ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim1.fastq.gz $WDIR/$OUT_flash/out.extendedFrags.fastq.gz
  elif [[ "$DATASET_libraryMethod" == UMI5RACEASYM ]]; then
    echo "Asymmetric sequencing dataset: looking for primers in Read1."
    cutadapt $args --trim-n -o $DATANAME.trim1.fastq.gz $WDIR/$INDIR/$DATA1
  elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(variable|variableNano)$ ]]; then
    echo "NEB-adaptored 5'RACE library with UMI's ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim1.fastq.gz $WDIR/$OUT_flash/out.extendedFrags.fastq.gz
  elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "Hinge dataset from NEB-adaptored library with UMI's (stitched reads from long amplicons) ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim1.fastq.gz $WDIR/$OUT_flash/out.extendedFrags.fastq.gz
  elif [[ "$DATASET_libraryType" == HINGE ]]; then
    echo "Hinge dataset: skipping 5' primer trimming ..."
    cp $WDIR/$OUT_flash/out.extendedFrags.fastq.gz $DATANAME.trim1.fastq.gz
  else
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --untrimmed-output=$DATANAME.untrimmed1.fastq.gz --trim-n -o $DATANAME.trim1.fastq.gz $WDIR/$OUT_flash/out.extendedFrags.fastq.gz
  fi
  unset args

  echo "####################################"
  time_msg "Running cutadapt for 3' end..."
  args=$(tr "\n" " " <$WDIR/$SCRDIR/adapter3.conf)
  if [[ "$DATASET_libraryMethod" == multiplexNEB ]] && [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "Hinge dataset from NEB-adaptored library: trimming the reverse primer sequences ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim2.fastq.gz $DATANAME.trim1.fastq.gz
  elif [[ "$DATASET_libraryMethod" == UMI5RACEASYM ]]; then
    echo "Asymmetric sequencing dataset: looking for primers in Read2."
    cutadapt $args --trim-n -o $DATANAME.trim2.fastq.gz $WDIR/$INDIR/$DATA2
  elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(variable|variableNano)$ ]]; then
    echo "NEB-adaptored 5'RACE library with UMI's ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim2.fastq.gz $DATANAME.trim1.fastq.gz
  elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "Hinge dataset from NEB-adaptored library with UMI's (stitched reads from long amplicons) ..."
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --trim-n -o $DATANAME.trim2.fastq.gz $DATANAME.trim1.fastq.gz
  else
    cutadapt $args -m $MINLENGTH -M $MAXLENGTH --untrimmed-output=$DATANAME.untrimmed2.fastq.gz --trim-n -o $DATANAME.trim2.fastq.gz $DATANAME.trim1.fastq.gz
  fi
  unset args

  ## accounting
  buildReadAccountingSummary cutadaptStepAcct
  buildUniqueAccountingSummary cutadaptStepAcct
}

# function: the fastx step
# arguments: libraryMethod revprimer libraryType
function fastxStep (){
  libraryMethod=$1
  revprimer=$2
  libraryType=$3

  fwdprimer='UMI5RACE'
  preamble='[ATCGN]{0,4}'
  barcode=${UMIbarcode:?} #'T[ATCG]{4}T[ATCG]{4}T[ATCG]{4}T'
  post='CTTG{1,7}'

  cd $WDIR/$OUT_fastxtk || { error "Error: FASTx output directory not accessible!"; }
  echo "#################################################"
  time_msg "Performing dataset collapsing steps ..."

  if [[ "$libraryMethod" == UMI5RACE ]]; then
    echo "Working with a $libraryMethod $libraryType library with directional adaptoring ..."
    echo "Processing UMI barcodes ..."
    ${zcat:?} $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz | fastq_to_fasta -Q 33 -v -n -o $DATANAME.trimmed.fasta
    perl $WDIR/$SCRDIR/fasta_barcode_count.pl $DATANAME.trimmed.fasta $preamble $barcode $post> $DATANAME.trimmed.bc_annot.fasta
    python3 $WDIR/$SCRDIR/fasta_barcode_consensus.py $DATANAME.trimmed.bc_annot.fasta > $DATANAME.trimmed.bc_annot.consensus.fastq
    fastq_to_fasta -Q 33 -v -n -i $DATANAME.trimmed.bc_annot.consensus.fastq -o $DATANAME.trimmed.bc_annot.consensus.fasta
    time_msg "Consensus building collapsed the set to" "`${grep:?} -c ">" $DATANAME.trimmed.bc_annot.consensus.fasta`" "sequences."
    echo "Unrecognized barcodes found in" "`$grep -c "barcode=unknown" $DATANAME.trimmed.bc_annot.fasta`" "sequences."
    echo "Cleaning up the basecalls ..."
    fastx_clipper -v -a N -i $DATANAME.trimmed.bc_annot.consensus.fasta -o $DATANAME.trimmed.bc_annot.consensus.noN.fasta
    ## clean up the input for igblast while retaining the UMI collapse
    cat $DATANAME.trimmed.bc_annot.consensus.noN.fasta| sed "s|;barcode=[^[:space:]]*;retained=|\-|" > $DATANAME.trimmed.bc_annot.consensus.noN.clean.fasta
    cp $DATANAME.trimmed.bc_annot.consensus.noN.clean.fasta $WDIR/$OUT_igblast/input.fasta

  elif [[ "$libraryMethod" == UMI5RACEASYM ]]; then
    echo "Working with a $libraryMethod $libraryType library adaptored with the NEB kit ..."
    echo "Recognizing/transferring UMI barcodes from asymmetric sequencing reads ..."
    echo "Forward primer name: \"$fwdprimer\"; reverse primer name: \"$revprimer\""
    $zcat $WDIR/$OUT_cutadapt/$DATANAME.trim1.fastq.gz > $DATANAME.trim1.fastq
    $zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz > $DATANAME.trim2.fastq
    perl $WDIR/$SCRDIR/fastq_asym_barcode_transfer.pl $DATANAME.trim1.fastq $DATANAME.trim2.fastq $fwdprimer $revprimer "$preamble" "$barcode" "$post" > $DATANAME.trim1.bc_annot.fastq
    rm $DATANAME.trim1.fastq $DATANAME.trim2.fastq
    echo "Ordering reads using UMI barcodes ..."
    perl $WDIR/$SCRDIR/fastq_asym_barcode_order.pl $DATANAME.trim1.bc_annot.fastq > $DATANAME.trim1.bc_annot.ordered.fastq
    echo "Trimming reads to quality of 15."
    cutadapt -q 15 -o $DATANAME.trim1.bc_annot.ordered_q15.fastq $DATANAME.trim1.bc_annot.ordered.fastq
    echo "Calculating consensus sequences for UMI-barcoded read clusters ..."
    python3 $WDIR/$SCRDIR/fastq_barcode_consensus.py $DATANAME.trim1.bc_annot.ordered_q15.fastq --min_size 2 > $DATANAME.trim1.bc_annot.ordered.cons.fastq
    time_msg "Consensus building collapsed the set to" "`${grep:?} -c "^@MIG" $DATANAME.trim1.bc_annot.ordered.cons.fastq`" "sequences."
    echo "Unrecognized barcodes found in" "`$grep -c "barcode=unknown" $DATANAME.trim1.bc_annot.fastq`" "sequences."
    perl $WDIR/$SCRDIR/fastq_barcode_consensus_interleaved_filter.pl $DATANAME.trim1.bc_annot.ordered.cons.fastq > $DATANAME.trim1.bc_annot.ordered.cons.interleaved.fastq

    echo "Performing FLASH to rebuild the amplicons from UMI cluster consensus sequences."
    asymmetricSequencingExtension "400" ${FLASH_minoverlap:?} "0.5" "$WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.cons.interleaved.fastq"
    $zcat $WDIR/$OUT_flash/UMI5RACEASYM.extendedFrags.fastq.gz| fastq_to_fasta -Q 33 -v -o $DATANAME.UMIcluster.extended.fasta
    touch $DATANAME.UMIcluster.extended.fasta # in case the previous command produced an empty file
    cp $DATANAME.UMIcluster.extended.fasta $WDIR/$OUT_igblast/input.fasta

  elif [[ "$libraryMethod" == UMI5RACENEB ]] && [[ "$libraryType" =~ ^(HINGE|HINGENano|variable|variableNano)$ ]]; then
    echo "Working with a $libraryMethod $libraryType library adaptored with the NEB kit ..."
    echo "Fixing orientation of the reverse reads in the symmetrically adaptored dataset..."
    echo "Forward primer name: \"$fwdprimer\"; reverse primer name: \"$revprimer\""
    $zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz | fastq_to_fasta -Q 33 -v -n -o $DATANAME.trimmed.fasta
    perl $WDIR/$SCRDIR/fastx_asym_orientation_fix.pl $DATANAME.trimmed.fasta $fwdprimer $revprimer >$DATANAME.trimmed.allorient.fasta
    echo "Dropping the amplicons with unknown orientation..."
    $grep -v ";orient_unk" $DATANAME.trimmed.allorient.fasta| $grep -A1 ">" | $grep -v "\-\-" > $DATANAME.trimmed.orient.fasta
    echo "Uknown orientations (primer recognition problems) encountered: `$grep -c ";orient_unk" $DATANAME.trimmed.allorient.fasta` times"
    echo "Ordering reads using UMI barcodes from FLASH-extended data..."
    perl $WDIR/$SCRDIR/fasta_barcode_count.pl $DATANAME.trimmed.orient.fasta "$preamble" "$barcode" "$post"> $DATANAME.trimmed.orient.bc_annot.fasta

    if [[ "$libraryType" =~ ^(variableNano|HINGE|HINGENano)$ ]]; then
       # This sequence is expected to be stitched; remove the forward read (before the stitch)
       echo "Releasing the 3' sequence from the stitched reads..."
       # cutadapt: the '-N' switch turns off the wildcard matching in the adapter sequences
       #           since the stitch contains a string of N's, we need to use this.
       cutadapt -O 25 -e 0 -N -g $READ_stitch $DATANAME.trimmed.orient.bc_annot.fasta -o $DATANAME.trimmed.orient.bc_annot.3prime.fasta
       # needed for proper MIG analysis
       cp $DATANAME.trimmed.orient.bc_annot.3prime.fasta $DATANAME.trimmed.orient.bc_annot.ordered.fasta

        echo "Determine the consensus sequence..."
       python3 $WDIR/$SCRDIR/fasta_barcode_consensus.py $DATANAME.trimmed.orient.bc_annot.3prime.fasta > $DATANAME.trimmed.orient.bc_annot.3prime.consensus.fastq
       fastq_to_fasta -Q 33 -v -n -i $DATANAME.trimmed.orient.bc_annot.3prime.consensus.fastq -o $DATANAME.trimmed.orient.bc_annot.3prime.consensus.fasta
       time_msg "Consensus building collapsed the set to" "`$grep -c ">" $DATANAME.trimmed.orient.bc_annot.3prime.consensus.fasta`" "sequences."
       echo "Unrecognized barcodes found in" "`$grep -c "barcode=unknown" $DATANAME.trimmed.orient.bc_annot.3prime.fasta`" "sequences."
       cp $DATANAME.trimmed.orient.bc_annot.3prime.consensus.fasta $WDIR/$OUT_igblast/input.fasta
    else
       # This sequence should be properly extended.
       echo "Determine the consensus sequence..."
       python3 $WDIR/$SCRDIR/fasta_barcode_consensus.py $DATANAME.trimmed.orient.bc_annot.fasta > $DATANAME.trimmed.orient.bc_annot.consensus.fastq
       fastq_to_fasta -Q 33 -v -n -i $DATANAME.trimmed.orient.bc_annot.consensus.fastq -o $DATANAME.trimmed.orient.bc_annot.consensus.fasta
       time_msg "Consensus building collapsed the set to" "`$grep -c ">" $DATANAME.trimmed.orient.bc_annot.consensus.fasta`" "sequences."
       echo "Unrecognized barcodes found in" "`$grep -c "barcode=unknown" $DATANAME.trimmed.orient.bc_annot.fasta`" "sequences."
       # needed for proper MIG accounting
       cp $DATANAME.trimmed.orient.bc_annot.fasta $DATANAME.trimmed.orient.bc_annot.ordered.fasta
       cp $DATANAME.trimmed.orient.bc_annot.consensus.fasta $WDIR/$OUT_igblast/input.fasta
    fi

  elif [[ "$libraryMethod" == multiplexNEB ]] && [[ "$libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "Working with a $libraryMethod $libraryType library adaptored with the NEB kit ..."
    fwdprimer='multiplex'
    echo "Fixing orientation of the reverse reads in the symmetrically adaptored HINGE dataset..."
    echo "Forward primer name: \"$fwdprimer\"; reverse primer name: \"$revprimer\""
    $zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz | fastq_to_fasta -Q 33 -v -n -o $DATANAME.trimmed.noN.fasta
    perl $WDIR/$SCRDIR/fastx_asym_orientation_fix.pl $DATANAME.trimmed.noN.fasta $fwdprimer $revprimer >$DATANAME.trimmed.noN.orient.fasta
    echo "Dropping the amplicons with unknown orientation..."
    $grep -v ";orient_unk" $DATANAME.trimmed.noN.orient.fasta| $grep -A1 ">" | $grep -v "\-\-" > $DATANAME.trimmed.noN.orient.clean.fasta
    echo "Uknown orientations (primer recognition problems) encountered: `$grep -c ";orient_unk"  $DATANAME.trimmed.noN.orient.fasta` times"
    time_msg "Running fastx_collapser ..."
    fastx_collapser -Q 33 -v -i $DATANAME.trimmed.noN.orient.clean.fasta -o $DATANAME.trimmed.noN.collapsed.fasta
    cp $DATANAME.trimmed.noN.collapsed.fasta $WDIR/$OUT_igblast/input.fasta

  else
    error "This pipeline doesn't know how to treat this library: \n$libraryType prepared using $libraryMethod"
  fi
  ## accounting
  buildReadAccountingSummary fastxStepAcct
  buildUniqueAccountingSummary fastxStepAcct
}

# function: the igblast step
# arguments: species
function IgBLASTstep (){
  species=$1

  case $species in
    Hs)
      echo "Setup: Human dataset $DATANAME."
      IGBLAST_species="human"
      ;;
    Rh)
      echo "Setup: Macaque dataset $DATANAME."
      IGBLAST_species="rhesus_monkey"
      ;;
    Mm)
      echo "Setup: Murine dataset $DATANAME."
      IGBLAST_species="mouse"
      ;;
    *)
      error "$species species not defined! Check that the dataset name follows correct nomenclature."
  esac

  cd $WDIR/$OUT_igblast || { error "Error: IgBLAST output directory not accessible!"; }

  echo "###############################"
  time_msg "Running igblastn ..."

  # UMI's are too far upstream of constant-region target sequences to yield enough
  #    variable-region sequence for a meaningful IgBLAST run ... skip this
  if [[ "$libraryMethod" == UMI5RACENEB ]] && [[ "$libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "UMI-tagged HINGE amplicons should not have variable-region sequences ... Skipping the igblastn step."
    touch $DATANAME.aa.igblast_out
    return 0
  fi

  echo "IgBLAST options: $species sequences ... Using $IGBLAST_species database."
  echo "Sample commandline: "
  echo "### igblastn -organism $IGBLAST_species"
  echo "###          -germline_db_V $IGDATA/database/${IGBLAST_species}_gl_V"
  echo "###          -germline_db_D $IGDATA/database/${IGBLAST_species}_gl_D"
  echo "###          -germline_db_J $IGDATA/database/${IGBLAST_species}_gl_J"
  echo "###          -auxiliary_data $IGDATA/optional_file/${IGBLAST_species}_gl.aux"
  echo "###          -show_translation"
  echo "###          -query input.fasta"
  echo "###          -num_threads ${IGBLAST_numthreads:?}"
  echo "###          -out $DATANAME.igblast_out"

  echo "Splitting the input file into 100,000-sequence blocks."
  split --verbose --lines=200000 input.fasta input_fasta_split.

  igblast_STARTTIME=$(date +%s)

  for f in input_fasta_split.*; do
    g=${f#*.}
    igblastn -organism $IGBLAST_species \
             -germline_db_V $IGDATA/database/${IGBLAST_species}_gl_V \
             -germline_db_D $IGDATA/database/${IGBLAST_species}_gl_D \
             -germline_db_J $IGDATA/database/${IGBLAST_species}_gl_J \
             -auxiliary_data $IGDATA/optional_file/${IGBLAST_species}_gl.aux \
             -show_translation \
             -query $f \
             -num_threads $IGBLAST_numthreads -out $DATANAME.${g}.igblast_out
    time_msg "Completed IgBLAST annotation of $f"
  done

  igblast_ENDTIME=$(date +%s)
  echo "The igblast step took $[$igblast_ENDTIME - $igblast_STARTTIME] seconds to complete."
}

# function: the igblast output processing
# arguments: chain, primer, libraryMethod, libraryType
function IgBLASToutputProcessing (){
  chain=$1
  primer=$2
  libraryMethod=$3
  libraryType=$4

  cd $WDIR/$OUT_igblast || { error "Error: IgBLAST output directory not accessible!"; }

  echo "######################################"
  time_msg "Processing igblastn output..."

  # UMI's are too far upstream of constant-region target sequences to yield enough
  #    variable-region sequence for a meaningful IgBLAST run ... skip this
  if [[ "$libraryMethod" == UMI5RACENEB ]] && [[ "$libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
    echo "We skipped the igblastn step for this $libraryType data ... Skipping the output processing..."
    cp input.fasta $DATANAME.igblast.prod.scrub.clon.fasta
    return 0
  fi

  for f in input_fasta_split.*; do
    g=${f#*.}
    if [[ -f $DATANAME.${g}.igblast_out ]]; then
      python3 $WDIR/$SCRDIR/igblast-out_harvester.py $DATANAME.${g}.igblast_out $f \
      >> $DATANAME.igblast.fasta
      time_msg "Completed transferring annotations from $DATANAME.${g}.igblast_out"
      rm $f # clean up the split-up fasta files
    else
      echo "Error!!! The file $DATANAME.${g}.igblast_out is missing. IgBLAST annotation was not completed."
    fi
  done

  # remove improperly truncated sequences and those containing stop codons in the CDR3aa (there may still be stops in the rest of the sequence!!!)
  if [[ "$libraryType" =~ ^(HINGE|HINGENano)$ ]]; then
      echo "Sorting out the productively-rearranged sequences for a hinge dataset (expecting truncations) ..."
      $grep -v "CDR3aa:\w*\*\w*" $DATANAME.igblast.fasta| $grep -A 1 "In-frame\tYes\t[\+-]\t[^\t]*\t?Vtruncated" | $grep -v "\-\-" > $DATANAME.igblast.prod.fasta
  elif [[ "$libraryType" == variableNano ]]; then
      echo "Sorting out the productively-rearranged sequences for a variableNano dataset (expecting truncations) ..."
      $grep -v "CDR3aa:\w*\*\w*" $DATANAME.igblast.fasta| $grep -A 1 "In-frame\tYes\t[\+-]\t[^\t]*\t?Vtruncated" | $grep -v "\-\-" > $DATANAME.igblast.prod.fasta
  else
      echo "Sorting out the productively-rearranged sequences..."
      $grep -v "CDR3aa:\w*\*\w*" $DATANAME.igblast.fasta| $grep -A 1 "Yes\t[\+-]\t[^\t]*\t?Vintact" | $grep -v "\-\-" > $DATANAME.igblast.prod.fasta
  fi

  echo "Removing the invalid (or unrecognized) sequences..."
  if [[ "$chain" =~ ^(IgM|IgG)$ ]];  then
    $grep -v "CDR3aa:0null0" $DATANAME.igblast.prod.fasta | $grep -A 1 "\t\VH\t" | $grep -v "\-\-" > $DATANAME.igblast.prod.scrub.fasta
  elif [[ "$chain" == IgK ]]
  then
    $grep -v "CDR3aa:0null0" $DATANAME.igblast.prod.fasta | $grep -A 1 "\t\VK\t" | $grep -v "\-\-" > $DATANAME.igblast.prod.scrub.fasta
  elif [[ "$chain" == IgL ]]
  then
    $grep -v "CDR3aa:0null0" $DATANAME.igblast.prod.fasta | $grep -A 1 "\t\VL\t" | $grep -v "\-\-" > $DATANAME.igblast.prod.scrub.fasta
  fi

  echo "Generating a clonotype dictionary..."
  perl $WDIR/$SCRDIR/clonotype_count.pl $DATANAME.igblast.prod.scrub.fasta > $DATANAME.igblast.prod.scrub.clonotype_dict
  $grep -v Vambig $DATANAME.igblast.prod.scrub.clonotype_dict |$grep -v "\w{5}\t\w+\t\d+\t[1234]\t" > $DATANAME.igblast.prod.scrub.5up.clonotype_dict

  echo "Annotating the productively-rearranged sequences with predicted clonotype information..."
  perl $WDIR/$SCRDIR/clonotype_annotate.pl $DATANAME.igblast.prod.scrub.clonotype_dict $DATANAME.igblast.prod.scrub.fasta > $DATANAME.igblast.prod.scrub.clon.fasta

  echo "Counting sequences..."
  $grep -c ">" *.fasta

  if [[ "$libraryMethod" =~ ^(UMI5RACEASYM|UMI5RACENEB)$ ]]; then
    sed -E '/^>/s/MIG([[:digit:]]+)[^[:space:]]+;retained=/\1\-/g' $DATANAME.igblast.prod.scrub.clon.fasta > $DATANAME.fastaviewer.fasta
  else
    cp $DATANAME.igblast.prod.scrub.clon.fasta $DATANAME.fastaviewer.fasta
  fi

  head -2 $DATANAME.igblast.prod.scrub.clonotype_dict| sed "s|^# ||"
  echo "Chimera tagged: " "`$grep -c "\-chimera" $DATANAME.igblast.prod.scrub.clon.fasta`"

  echo "Unrecognized CDR3's in the productively-rearranged set: " "`$grep "CDR3aa:0null0" -c $DATANAME.igblast.prod.fasta`"

  ## accounting
  buildReadAccountingSummary IgBLASTstepAcct
  buildUniqueAccountingSummary IgBLASTstepAcct
}

# function: the hinge processing step
# arguments: species
function hingeProcessingStep (){
  species=$1

  case $species in
    Hs)
      echo "Setup: Human dataset $DATANAME."
      BLAST_species="human"
      ;;
    Rh)
      echo "Setup: Macaque dataset $DATANAME."
      BLAST_species="rhesus_monkey"
      ;;
    *)
      error "$species species not defined! Check that the dataset name follows correct nomenclature."
  esac

  cd $WDIR/$OUT_igblast || { error "Error: IgBLAST output directory not accessible!"; }

  echo "########################################"
  time_msg "Processing hinge data ..."
  echo "Assigning subtype using blastn ..."
  echo "BLAST options: $species sequences ... Using $BLAST_species database."
  echo "commandline: "
  echo "### blastn -db $BLAST_DATA/${BLAST_species}_constant "
  echo "###        -query $DATANAME.igblast.prod.scrub.clon.fasta"
  echo "###        -num_threads ${BLAST_numthreads:?}"
  echo "###        -outfmt \"7 qseqid sseqid pident qlen slen length qcovs bitscore evalue\""
  echo "###        -max_target_seqs 3"
  echo "###        -out $DATANAME.igblast.prod.scrub.clon.blast_out"

  echo "Splitting the $DATANAME.igblast.prod.scrub.clon.fasta file into 100,000-sequence blocks."
  split --verbose --lines=200000 ${DATANAME}.igblast.prod.scrub.clon.fasta \
            ${DATANAME}_igblast_prod_scrub_clon_fasta_split.

  blastn_STARTTIME=$(date +%s)
  for f in ${DATANAME}_igblast_prod_scrub_clon_fasta_split.*; do
    g=${f#*.}
    blastn -db $BLAST_DATA/${BLAST_species}_constant \
           -query $f \
           -num_threads $BLAST_numthreads \
           -outfmt "7 qseqid sseqid pident qlen slen length qcovs bitscore evalue" \
           -max_target_seqs 3 \
           -out $DATANAME.igblast.prod.scrub.clon.${g}.blast_out
    time_msg "Completed HINGE BLAST annotation of $f"
  done
  blastn_ENDTIME=$(date +%s)
  echo "The hinge blast step took $[$blastn_ENDTIME - $blastn_STARTTIME] seconds to complete."

  echo "########################################"
  time_msg "Parsing the hinge blast output..."

  for f in ${DATANAME}_igblast_prod_scrub_clon_fasta_split.*; do
    g=${f#*.}
    if [[ -f $DATANAME.igblast.prod.scrub.clon.${g}.blast_out ]]; then
      perl $WDIR/$SCRDIR/hinge_blast_out_harvester.pl \
           $DATANAME.igblast.prod.scrub.clon.${g}.blast_out $f \
      >> $DATANAME.igblast.prod.scrub.clon.subclass.fasta
      time_msg "Completed transferring annotations from $DATANAME.${g}.igblast_out"
      rm $f # clean up the split-up fasta files
    else
      echo "Error!!! The file $DATANAME.${g}.igblast_out is missing. IgBLAST annotation was not completed."
    fi
  done

  echo "Generating a representative subset for each clonal cluster..."
  perl $WDIR/$SCRDIR/subclass_subset.pl $DATANAME.igblast.prod.scrub.clon.subclass.fasta > $DATANAME.igblast.prod.scrub.clon.subclass.subset.fasta

  echo "Preparing a CDR3aa dictionary..."
  $grep ">" $DATANAME.igblast.prod.scrub.clon.subclass.subset.fasta | sed "s|^>[[:digit:]]\+\-[[:digit:]]\+;\([[:alnum:]]\+\-[[:digit:]]\+\)\S*\t.\+CDR3aa\:\([[:alpha:]]\+\)\t.\+|>\1\t\2|" |sort |uniq |sed "s|\t|\n|" > $DATANAME.igblast.prod.scrub.clon.subclass.subset.CDR3aa_dict.fasta

}

# function: accounting of reads passing through the steps of the pipeline
# arguments: accounting step
function buildReadAccountingSummary (){
  step=$1

  case $step in
    initialize )
      echo "# Pipeline progress accounting for processed reads:" > $WDIR/$OUTDIR/$DATANAME.accounting.csv
      echo "\"Total\"," "`$zcat $WDIR/$INDIR/$DATA1| $grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      ;;
    extendStepAcct )
      echo "\"Extended\"," "`$zcat $WDIR/$OUT_flash/out.extendedFrags.fastq.gz| $grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      ;;
    cutadaptStepAcct )
      if [[ "$DATASET_libraryMethod" =~ ^(UMI5RACEASYM|UMI5RACENEB|multiplexNEB)$ ]]; then
        echo "\"5trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim1.fastq.gz| $grep -v "^@.+;no_adapter$" |$grep -c "^@.+;.+_adapter$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"3trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz| $grep -v "^@.+;no_adapter$" |$grep -c "^@.+;.+_adapter$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      else
        echo "\"5trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim1.fastq.gz|$grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"3trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz|$grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      fi
      ;;
    fastxStepAcct )
      if [[ "$DATASET_libraryMethod" == UMI5RACE ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.bc_annot.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"clean\"," "`$grep "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.bc_annot.consensus.noN.fasta | cut -d "=" -f4| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(variableNano|HINGE|HINGENano)$ ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.3prime.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"clean\"," "`$grep "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.3prime.consensus.fasta | cut -d "=" -f4| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv

      elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" == variable ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"clean\"," "`$grep "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.consensus.fasta | cut -d "=" -f4| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv

      elif [[ "$DATASET_libraryMethod" == UMI5RACEASYM ]]; then
        echo "\"barcoded\"," "`$grep -v "^@.+;barcode=unknown$" $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.fastq| $grep -c "^@.+;barcode="`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"valid\"," "`$grep -c ";valid;" $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.fastq`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"paired\"," "`$grep "^@.+barcode=" $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.cons.interleaved.fastq| uniq | cut -d ";" -f4| cut -d "=" -f2| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"AsymmetricExt\"," "`$zcat $WDIR/$OUT_flash/UMI5RACEASYM.extendedFrags.fastq.gz|$grep "^@.+barcode="| cut -d ";" -f4| cut -d "=" -f2| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"clean\"," "`$grep "^>" $WDIR/$OUT_fastxtk/$DATANAME.UMIcluster.extended.fasta | cut -f1| cut -d "=" -f3| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      else
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.noN.fasta`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      fi
      ;;
    IgBLASTstepAcct )
      if [[ "$DATASET_libraryMethod" =~ ^(UMI5RACEASYM|UMI5RACENEB)$ ]] && [[ "$DATASET_libraryType" =~ ^(variable|variableNano)$ ]]; then
        echo "\"productive\"," "`$grep "^>" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta | cut -f1| cut -d "=" -f3| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"Chimera\"," "`$grep "^>.+\-chimera" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta| cut -f1| cut -d "=" -f3| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      else
        echo "\"productive\"," "`$grep "^>" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta | cut -d ";" -f1| cut -d "-" -f2| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
        echo "\"Chimera\"," "`$grep "^>.+\-chimera" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta| cut -d ";" -f1| cut -d "-" -f2| awk 'BEGIN{s=0}{s+=$1}END{print s}'`" >> $WDIR/$OUTDIR/$DATANAME.accounting.csv
      fi
  esac
}

# function: accounting of unique sequences passing through the steps of the pipeline
# arguments: accounting step
function buildUniqueAccountingSummary (){
  step=$1

  case $step in
    initialize )
      echo "# Pipeline progress accounting for processed unique sequences:" > $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      echo "\"Total\"," "`$zcat $WDIR/$INDIR/$DATA1| $grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      ;;
    extendStepAcct )
      echo "\"Extended\"," "`$zcat $WDIR/$OUT_flash/out.extendedFrags.fastq.gz| $grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      ;;
    cutadaptStepAcct )
      if [[ "$DATASET_libraryMethod" =~ ^(UMI5RACEASYM|UMI5RACENEB|multiplexNEB)$ ]]; then
        echo "\"5trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim1.fastq.gz| $grep -v "^@.+;no_adapter$" |$grep -c "^@.+;.+_adapter$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"3trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz| $grep -v "^@.+;no_adapter$" |$grep -c "^@.+;.+_adapter$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      else
        echo "\"5trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim1.fastq.gz|$grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"3trim\"," "`$zcat $WDIR/$OUT_cutadapt/$DATANAME.trim2.fastq.gz|$grep -c "^\+$"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      fi
      ;;
    fastxStepAcct )
      if [[ "$DATASET_libraryMethod" == UMI5RACE ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.bc_annot.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.bc_annot.consensus.noN.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" =~ ^(HINGE|HINGENano|variableNano)$ ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.3prime.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.3prime.consensus.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      elif [[ "$DATASET_libraryMethod" == UMI5RACENEB ]] && [[ "$DATASET_libraryType" == variable ]]; then
        echo "\"barcoded\"," "`$grep -v "barcode=unknown" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.fasta| $grep -c "^>"`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.consensus.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      elif [[ "$DATASET_libraryMethod" == UMI5RACEASYM ]]; then
        echo "\"barcoded\"," "`$grep -v "^@.+;barcode=unknown$" $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.fastq| $grep -c "^@.+;barcode="`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"paired MIGs\"," "`$grep -c "^@.+barcode=" $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.cons.interleaved.fastq`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"AsymmetricExt\"," "`$zcat $WDIR/$OUT_flash/UMI5RACEASYM.extendedFrags.fastq.gz|$grep -c "^@.+barcode="`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.UMIcluster.extended.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      else
        echo "\"clean\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.noN.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"collapsed\"," "`$grep -c "^>" $WDIR/$OUT_fastxtk/$DATANAME.trimmed.noN.collapsed.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      fi
      ;;
    IgBLASTstepAcct )
      if [[ "$DATASET_libraryMethod" == UMI5RACEASYM ]]; then
        echo "\"productive\"," "`$grep -c "^>" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"Chimera\"," "`$grep -c "^>.+\-chimera" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      else
        echo "\"productive\"," "`$grep -c "^>" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
        echo "\"Chimera\"," "`$grep -c "^>.+\-chimera" $WDIR/$OUT_igblast/$DATANAME.igblast.prod.scrub.clon.fasta`" >> $WDIR/$OUTDIR/$DATANAME.uniqaccounting.csv
      fi
  esac
}
