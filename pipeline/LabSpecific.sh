#!/bin/bash
# collection of valid values relevant for Sather Lab
# shellcheck disable=SC2034 # Unused variable are invoked in other scripts

declare -A validSpecies=( ["Hs"]=1 ["Rh"]=1 ["Mm"]=1 )
declare -A validChains=( ["IgG"]=1 ["IgK"]=1 ["IgL"]=1 ["IgM"]=1 )
declare -A validAdjuvants=( ["MPLA_QS21"]=1 )
declare -A validDiseases=( ["HIV"]=1 ["malaria"]=1 )
declare -A validImmunogens=( ["PyCSP_NxA"]=1 )
declare -A validPrimers=( ["STD"]=1 ["HsRhIgGCH1Rev"]=1 ["HuIgGint"]=1 ["vv874pool876"]=1 )
declare -A validResponsiblePersons=( ["NSather"]=1 ["VVigdorovich"]=1 )
declare -A validTimepoints=(["PB"]=1 ["m2wk"]=1 ["p0wk"]=1 ["p2wk"]=1 ["p4wk"]=1 ["p6wk"]=1 )
declare -A validTrials=( ["IMRAS"]=1 )
declare -A validTissues=(["spleen"]=1 ["marrow"]=1 ["blood"]=1 ["node"]=1 ["PBMC"]=1 )

# primers [species:chain:libraryMethod] leading to the file located in adapters/*
declare -A primers5=(["all:all:UMI5RACE"]="adapter5UMI.conf"
  ["all:all:UMI5RACEASYM"]="adapterUMIasym.conf"
  ["all:all:UMI5RACENEB"]="adapter5UMIasym.conf"
  ["all:all:multiplexNEB"]="adapter5multiplexNEB.conf"
)

declare -A primers3=(["all:all:UMI5RACEASYM"]="adapterUMIasym.conf"
  ["all:all:UMI5RACENEB"]="adapter3UMIasym.conf"
  ["all:all:multiplexNEB"]="adapter3multiplexNEB.conf"
)

####################################################################
#  Edit the lab-specific custom settings above this line ...
####################################################################

# Process-affecting variables
# multiplexNEB -- multiplex 5' primer pool with Illumina adapters ligated on (NEB)
# UMI5RACE -- 5'RACE library using SMARTer adapter with UMIs,
#              Illumina adapters incorporated into primers
# UMI5RACEASYM -- 5'RACE library using SMARTer adapter with UMIs,
#                  Illumina adapters ligated on (NEB), asymmetric sequencing
# UMI5RACENEB -- 5'RACE library using SMARTer adapter with UMIs, Illumina adapters ligated
#                  on (NEB)

declare -A validLibraryMethods=( ["UMI5RACE"]=1 ["UMI5RACEASYM"]=1
  ["UMI5RACENEB"]=1 ["multiplexNEB"]=1
)

# variable -- variable-region libraries
# HINGE -- constant-region libraries
# *Nano -- "nanochip" reading out 150bp in each direction
#               (cheap, quick but not extendible by FLASH)

declare -A validLibraryTypes=(["variable"]=1 ["HINGE"]=1 ["variableNano"]=1
  ["HINGENano"]=1
)
