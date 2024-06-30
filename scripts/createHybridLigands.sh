#!/bin/bash
set -e

echo "CADD: Create Hybrid Ligands " >>$LOGFILE
# required external software
module use $CADD_SOFTWARE_MODULES
module load python/3.11.1
module load openbabel/3.1.1
module load ambertools/22
module load acpype/2022.7.21
module load gromacs/2023.2
module load rdkit/2023.3.1

# STaGE itself
module load stage/1.0.0
export GMXLIB=$GMXDATA

OUTPUT_PATH=${CADD_OUTPUT_DIR}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}
target=${CADD_TARGET}

cd ${OUTPUT_PATH}/${target}

echo "Creating hybrid ligands for ${target} ..."

# run ligand pairing script
cp ${PATH_TO_SCRIPTS}/CB_ligand_pairing_v4.py .
cp ${PATH_TO_SCRIPTS}/MCS_sebastian.py .

python3 CB_ligand_pairing_v4.py

if ! [ -f ligandPairs.json ]
then
echo "Ligand pairing was completely unsuccessful. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload openbabel/3.1.1
module unload ambertools/22
module unload acpype/2022.7.21
module unload gromacs/2023.2
module unload rdkit/2023.3.1

export GMXLIB=""

ligands1=($(grep "Ligand_1" ligandPairs.json | cut -d "\"" -f 4 | cut -d "/" -f 2))
ligands2=($(grep "Ligand_2" ligandPairs.json | cut -d "\"" -f 4 | cut -d "/" -f 2))

# clean up
rm -r CB_ligand_pairing_v4.py MCS_sebastian.py __pycache__

# merge topologies for pairs of ligands
for index in $(seq 0 $(( ${#ligands1[@]}-1 ))); do

mkdir edge_${ligands1[${index}]}_${ligands2[${index}]}
cd edge_${ligands1[${index}]}_${ligands2[${index}]}

for i in 0 1; do
for j in 0 1; do

mkdir pose_${i}_${j}
cd pose_${i}_${j}

sbatch --wait ${PATH_TO_SCRIPTS}/createHybridLigands2.sh &

cd ..
done
done

cd ..
done
wait
echo "CADD: Create Hybrid Ligands -- finished" >>$LOGFILE
