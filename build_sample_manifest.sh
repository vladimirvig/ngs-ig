#!/bin/bash
# Script to generate sample manifest from dataset path.
# expected order of fields (inherited from FASTAViewer SatherLabUtils.R parser:
# 1. Species (required first)
# 2. Subject/individual (required second)
# 3. Trial code
# 4. Immunogen
# 5. Adjuvant
# 6. Primers used
# 7. Timepoint
# 8. Disease (optional)
# 9. Tissue (optional)
# 10. Replicate (optional)
# 11. Responsible person (optional)
# 12. Library type (variable, HINGE)
# 13. (optional)
# 14. (optional)
# 15. (optional)
# 16. Library method (required)
# 17. Chain (required last)

#
# Expected name convention for processing:
# exp_date/species-...-chain
# everything in between is parsed against a database of valid values stored in LabSpecific.sh
#

if [[ $1 = '' ]]; then
	echo "commandline: build_sample_manifest.sh <directory_name>"
	exit 1
fi

absPath1=$1

if [[ -f $absPath1/scripts/LabSpecific.sh ]]; then
	source $absPath1/scripts/LabSpecific.sh
fi

# get experiment date
absPath2=`dirname $absPath1`
runDirectory=${absPath2##/*/}
experimentDate=${runDirectory%%-*}

# get dataset name
sampleName=`basename $absPath1`
sampleNameString=$sampleName # this string will be modified
spacer="+"

sampleSpecies=''
sampleChain=''
sampleLibraryMethod=''
sampleImmunogen=''
sampleAdjuvant=''
sampleDisease=''
samplePrimer=''
sampleResponsiblePerson=''
sampleTimepoint=''
sampleTrial=''
sampleTissue=''
sampleSubject=''
sampleReplicate=''
sampleLibraryType=''

# date format test 1: length
if [[ ${#experimentDate} != 8 ]]
then
  echo "Invalid date format for the experiment! (wrong length): $experimentDate"
  exit 1
fi
# date format test 2: numeric
if ! [[ $experimentDate =~ ^[0-9]+$ ]]
then
  echo "Invalid date format for the experiment! (non-numeric)"
  exit 1
fi

# harvest first (anchor) field: species
sampleSpecies=${sampleNameString%%-*}
sampleNameString=${sampleNameString#*-}
if ! [[ ${validSpecies[$sampleSpecies]} ]]
then
  echo "Invalid species detected in $sampleName: $sampleSpecies"
  exit 1
fi

# harvest last (anchor) field: chain
sampleChain=${sampleNameString##*-}
sampleNameString=${sampleNameString%-*}
if ! [[ ${validChains[$sampleChain]} ]]
then
  echo "Invalid chain detected in $sampleName: $sampleChain"
  exit 1
fi

# parse the rest of the fields that may have variable usage and order
IFS='-' read -r -a splitSampleNameArray <<< "$sampleNameString"
for i in "${!splitSampleNameArray[@]}"; do
  # check Library Method
  if ! [[ $sampleLibraryMethod ]] && [[ -n ${validLibraryMethods[${splitSampleNameArray[$i]}]} ]]; then
    sampleLibraryMethod=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Immunogen
  elif ! [[ $sampleImmunogen ]] && [[ ${validImmunogens[${splitSampleNameArray[$i]}]} ]]; then
    sampleImmunogen=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Adjuvant
  elif ! [[ $sampleAdjuvant ]] && [[ -n ${validAdjuvants[${splitSampleNameArray[$i]}]} ]]; then
    sampleAdjuvant=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Disease
  elif ! [[ $sampleDisease ]] && [[ -n ${validDiseases[${splitSampleNameArray[$i]}]} ]]; then
    sampleDisease=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Primer
  elif ! [[ $samplePrimer ]] && [[ -n ${validPrimers[${splitSampleNameArray[$i]}]} ]]; then
    samplePrimer=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Responsible Person
  elif ! [[ $sampleResponsiblePerson ]] && [[ -n ${validResponsiblePersons[${splitSampleNameArray[$i]}]} ]]; then
    sampleResponsiblePerson=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Timepoint
  elif ! [[ $sampleTimepoint ]] && [[ -n ${validTimepoints[${splitSampleNameArray[$i]}]} ]]; then
    sampleTimepoint=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
  # check Trial
  elif ! [[ $sampleTrial ]] && [[ -n ${validTrials[${splitSampleNameArray[$i]}]} ]]; then
    sampleTrial=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
	# check Tissue
  elif ! [[ $sampleTissue ]] && [[ -n ${validTissues[${splitSampleNameArray[$i]}]} ]]; then
    sampleTissue=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  #fi
	# check Library Type
	elif ! [[ $sampleLibraryType ]] && [[ -n ${validLibraryTypes[${splitSampleNameArray[$i]}]} ]]; then
    sampleLibraryType=${splitSampleNameArray[$i]}
    unset "splitSampleNameArray[$i]"
    i=0
  fi
done

if [[ ${#splitSampleNameArray[*]} -gt 1 ]]; then
  echo "# Warning: there are more fields in Sample Name than expected..."
  echo "# Warning: this usually means that unexpected values were encounered."
  echo "# I ended up with this after parsing: \"${splitSampleNameArray[*]}\""
  IFS=$spacer sampleSubject=${splitSampleNameArray[*]}
else
  sampleSubject=${splitSampleNameArray[*]}
fi

# this is a hack to keep old datasets cooperating
if ! [[ $sampleLibraryMethod ]]; then
	echo "# Warning: no libraryMethod resolved from sampleName."
	echo "# Warning: using multiplex as the default."
  sampleLibraryMethod='multiplex'
fi

##########
# Output
##########
echo "# Sample manifest for $sampleName from $experimentDate"
echo "processDate: `date +%Y%m%d_%H%M`"
echo "species: $sampleSpecies"
echo "subjectID: $sampleSubject"
echo "chain: $sampleChain"
echo "timepoint: $sampleTimepoint"
echo "immunogen: $sampleImmunogen"
echo "adjuvant: $sampleAdjuvant"
echo "primer: $samplePrimer"
echo "trial: $sampleTrial"
echo "replicate: $sampleReplicate"
echo "disease: $sampleDisease"
echo "tissue: $sampleTissue"
echo "libraryMethod: $sampleLibraryMethod"
echo "experimentDate: $experimentDate"
echo "libraryType: $sampleLibraryType"
echo "unused2: "
echo "unused3: "
echo "unused4: "
echo "responsiblePerson: $sampleResponsiblePerson"
