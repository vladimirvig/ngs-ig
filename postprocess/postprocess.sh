#!/bin/bash
#
#   postprocess.sh
#   Task: generate data visualization as a part of the pipeline
#  Warning: this is not a standalone script; it is meant to be called from
#      the main processing script because it depends on its environment setup
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

###################################################
echo "Summarizing run dataset retention into a graph ..."
cd $WDIR/$OUTDIR || { echo "Error: output directory not accessible!"; exit 1; }
Rscript $WDIR/$SCRDIR/postprocess/buildDataRetentionPlot.R --inputFile=$DATANAME.accounting.csv --outputFile=$DATANAME.retention.png
Rscript $WDIR/$SCRDIR/postprocess/buildDataRetentionPlot.R --inputFile=$DATANAME.uniqaccounting.csv --outputFile=$DATANAME.uniques.png

if [ $? -ne 0 ]; then
  echo "Warning: something went wrong with figure generation..."
fi

cd $WDIR/$OUTDIR || { echo "Error: output directory not accessible!"; exit 1; }

case ${DATASET_libraryMethod:?} in
  UMI5RACE)
    echo "Generating UMI group summary plots"
    cat $WDIR/${OUT_fastxtk:?}/$DATANAME.trimmed.bc_annot.fasta |${grep:?} -v "barcode=unknown" | $grep "^>.*;element=1$"| sed "s|^.*;size=||" | cut -d ";" -f1 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_retained_final.pdf
    cat $WDIR/$OUT_fastxtk/$DATANAME.trimmed.bc_annot.consensus.fasta |$grep -v "barcode=unknown" | $grep "^>MIG"| sed "s|^.*;retained=||" | cut -d ";" -f1 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_consensus_final.pdf
    ;;
  UMI5RACEASYM)
    echo "Generating UMI group summary plots..."
    echo "Examining IgBLAST output ..."
    cat $WDIR/${OUT_igblast:?}/$DATANAME.fastaviewer.fasta |$grep ">" | cut -d ";" -f1| cut -d "-" -f2 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_IgBLAST_final.pdf
    echo "Examining raw UMI groups ..."
    cat $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.fastq | $grep "^@.*;valid;.*;element=1$" | sed "s|^.*;size=||" | cut -d ";" -f1 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_retained_final.pdf
    echo "Examining interleaved consensus datasets..."
    cat $WDIR/$OUT_fastxtk/$DATANAME.trim1.bc_annot.ordered.cons.interleaved.fastq | $grep "^@MIG.*;retained=" | uniq | sed "s|^.*;retained=||" |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_paired_final.pdf
    ;;
  UMI5RACENEB)
    echo "Generating UMI group summary plots..."
    if [[ "${DATASET_libraryType:?}" =~ ^(variable|variableNano)$ ]]; then
      echo "Examining IgBLAST output ..."
      cat $WDIR/$OUT_igblast/$DATANAME.fastaviewer.fasta |$grep ">" | cut -d ";" -f1| cut -d "-" -f2 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
      mv output.pdf ${DATANAME}_IgBLAST_final.pdf
    fi
    echo "Examining raw UMI groups ..."
    cat $WDIR/$OUT_fastxtk/$DATANAME.trimmed.orient.bc_annot.ordered.fasta | $grep "^>.*;element=1$" | $grep -v ";barcode=unknown;" | sed "s|^.*;size=||" | cut -d ";" -f1 |Rscript $WDIR/$SCRDIR/postprocess/makeAbundanceFig.R
    mv output.pdf ${DATANAME}_retained_final.pdf
    ;;
  *)
    ;;
esac
