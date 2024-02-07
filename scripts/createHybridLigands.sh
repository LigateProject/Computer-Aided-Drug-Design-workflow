#!/bin/bash
set -e

# required external software
module load python/3.11.1
module load forSTaGE/openbabel/3.1.1
module load forSTaGE/ambertools/21
module load forSTaGE/acpype/2022.7.21
module load gromacs/2023.2 gromacs=gmx
module load rdkit/2023.3.1

# STaGE itself
module load stage/1.0.0
export GMXLIB=$GMXDATA

target=$(pwd | rev | cut -d "/" -f 1 | rev)

PATH_TO_SCRIPTS=CADD/scripts

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

ligands1=($(grep "Ligand_1" ligandPairs.json | cut -d "\"" -f 4 | cut -d "/" -f 2))
ligands2=($(grep "Ligand_2" ligandPairs.json | cut -d "\"" -f 4 | cut -d "/" -f 2))

# clean up
rm -r CB_ligand_pairing_v4.py MCS_sebastian.py __pycache__ ligandPairs.json

# merge topologies for pairs of ligands
for index in $(seq 0 $(( ${#ligands1[@]}-1 ))); do

mkdir edge_${ligands1[${index}]}_${ligands2[${index}]}
cd edge_${ligands1[${index}]}_${ligands2[${index}]}

#TODO: avoid repeating topology merge; can end up with different hybrid ligands between replicas now
for i in 0 1; do
for j in 0 1; do
# list input files
## pair pose 0 with pose 0 and pose 1 with pose 1
## but make sure that both ligands dominate the initial coordinates once for each pose pair
if (( j == 0 ))
then
file1=../../${ligands1[${index}]}/ligand_premature.itp
file2=../../${ligands2[${index}]}/ligand_premature.itp
file3=../../${ligands1[${index}]}/pose_${i}/ligand.mol2
file4=../../${ligands2[${index}]}/pose_${i}/ligand.mol2
file5=../../${ligands1[${index}]}/pose_${i}/ligand.gro
file6=../../${ligands2[${index}]}/pose_${i}/ligand.gro
elif (( j == 1 ))
then
file1=../../${ligands2[${index}]}/ligand_premature.itp
file2=../../${ligands1[${index}]}/ligand_premature.itp
file3=../../${ligands2[${index}]}/pose_${i}/ligand.mol2
file4=../../${ligands1[${index}]}/pose_${i}/ligand.mol2
file5=../../${ligands2[${index}]}/pose_${i}/ligand.gro
file6=../../${ligands1[${index}]}/pose_${i}/ligand.gro
fi

mkdir pose_${i}_${j}
cd pose_${i}_${j}

# merge topologies and .gro files
python3 $PATH_TO_SCRIPTS/mergeTopologies.py <<-eof
$file1
$file2
$file3
$file4
merged.itp
$file5
$file6
merged.gro
eof

# write topology summary
python3 $PATH_TO_SCRIPTS/writeTopologySummary.py

# fix structure of hybrid ligand
## currently, the hybrid ligand is ligand A in state 0 and ligand B in state 1
## therefore, the starting structure is the structure of ligand A with approriately placed dummy atoms needed for ligand B
python3 $PATH_TO_SCRIPTS/posResForLigandToFixStructure.py

mv merged.gro merged_old.gro
## make box 10 nm larger in each dimension for energy minimization
python3 $PATH_TO_SCRIPTS/changeLastLineOfGroFile.py <<-eof
merged_old.gro
10
eof
# The ligand in vacuum may have a charge, but potential electrostatic artifacts are irrelevant for this very short energy minimisation
gmx grompp -f $PATH_TO_SCRIPTS/mdp/fixStructureOfHybridLigand.mdp -c merged_old.gro -r merged_old.gro -p topol_ligandInWater.top -o merged.tpr -maxwarn 1
gmx mdrun -deffnm merged
rm merged_old.gro mdout.mdp merged.tpr merged.trr merged.edr merged.log

## go back to original box size
python3 $PATH_TO_SCRIPTS/changeLastLineOfGroFile.py <<-eof
merged.gro
-10
eof

# print .gro file for complex of protein and hybrid ligand
python3 $PATH_TO_SCRIPTS/printComplexGroFile.py <<-eof
../../conf.gro
merged.gro
full.gro
eof

# print position restraints file for ligand
python3 $PATH_TO_SCRIPTS/posResForLigand.py 

# include position restraints for protein correctly
for file in *; do
if grep -q "posre" ${file} && ! grep -q "posre_ligand" ${file}
then
sed -i ${file} -e "s=posre=../../posre=g"
fi
done

cd ..
done
done

cd ..
done

# clean up
rm topol.top conf.gro
for lig in $(ls -d */ | grep -v "edge"); do
rm -r ${lig}
done

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload forSTaGE/openbabel/3.1.1
module unload forSTaGE/ambertools/21
module unload forSTaGE/acpype/2022.7.21
module unload gromacs/2023.2
module unload rdkit/2023.3.1

export GMXLIB=""
