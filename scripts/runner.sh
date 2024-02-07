#!/bin/bash
set -e; shopt -s expand_aliases

script=runInputFileGeneration.sh

PATH_TO_SCRIPTS=$(pwd)
cd ${ligate}/CADDValidation
sbatch ${PATH_TO_SCRIPTS}/${script}
