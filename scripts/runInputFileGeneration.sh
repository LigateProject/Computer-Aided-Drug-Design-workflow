#!/bin/bash

#SBATCH --job-name=PSInputFileGeneration
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=8GB
#SBATCH -t 48:00:00
#SBATCH -o slurm-%J.out

set -e

OUTPUT_PATH=${ligate}/CADDValidation/
PATH_TO_SCRIPTS=CADD/scripts

targets=(bace_p2)

for target in ${targets[@]}; do

# work in /scratch not to overwhelm the file system
mkdir -p /scratch
cd /scratch
if [ -d ${target} ]
then
rm -rf ${target}
echo "The remainders of a previous unsuccessful workflow execution were removed."
fi
mkdir ${target}
cd ${target}

bash ${PATH_TO_SCRIPTS}/checkProtein.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/GROMACS_LiGen_integration.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/createHybridLigands.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/solvate.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/addIons.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/energyMinimisation.sh

# only deal with targets that were not removed by previous stages of the workflow
if ! [ -d ../${target} ]
then
continue
fi
bash ${PATH_TO_SCRIPTS}/prepareEquilibration.sh

# copy data to file system
cd ..
cp -r ${target} ${OUTPUT_PATH}
rm -r ${target}
echo "${target} successfully completed!"

done # targets

echo "DONE! :)"
