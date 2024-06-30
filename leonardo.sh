#!/bin/bash
#
#  ** CADD WORKFLOW LEONARDO **
#  Variables required for leonardo
#  source leonardo.sh before starting
#
#    Last update: 15/03/2024
#

# target list (separate targets by spaces)
export CADD_TARGET="AAT75224"

# number of ligands (TODO)
#export NUMBEROFLIGANDS=2

# change as appropriate
export CADD_HOME=/leonardo/home/userinternal/aemerson/ligate-workflows

# In case undefined
export SCRATCH=$CINECA_SCRATCH

# software environment  #

## software
### Change as appropriate but these directories are globally readable
#### local modules
export CADD_SOFTWARE=/leonardo/pub/userinternal/aemerson
export CADD_SOFTWARE_MODULES=${CADD_SOFTWARE}/modules


### CADD scripts - no  need to change
## where the scripts are - can be leonardo $HOME directory
export CADD_SCRIPTS_DIR=${CADD_HOME}/scripts

#--------------------------------------------------
#  run environment - change as necessary #

## where to run - should be $SCRATCH or $WORK (not $HOME)
export CADD_SUBMISSION_DIR=${PWD}
mkdir -p ${CADD_SUBMISSION_DIR}

# LIGEN results
export resultFile="${CADD_SUBMISSION_DIR}/${CADD_TARGET}/result.csv"
export mol2File="${CADD_SUBMISSION_DIR}/${CADD_TARGET}/docked.mol2"

# output - not in submission dir
export CADD_OUTPUT_DIR=${CADD_SUBMISSION_DIR}/Validation
mkdir -p ${CADD_OUTPUT_DIR}
export OUTPUT_PATH=${CADD_OUTPUT_DIR}


# inputs should be under $CADD_SUBMISSION_DIR
export CADD_INPUT_DIR=${CADD_SUBMISSION_DIR}

# workdir
export CADD_LOCAL_DIR=${CADD_SUBMISSION_DIR}/CADD_local
rm -rf ${CADD_LOCAL_DIR}
mkdir -p ${CADD_LOCAL_DIR}

# logging
export LOGFILE="$PWD/${CADD_TARGET}.log"

