#!/bin/bash

#SBATCH --job-name=createHybridLigand
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=32GB
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out
#SBATCH -A cin_staff

set -e

echo "CADD: Create Hybrid Ligands 2 ($(pwd))" >>$LOGFILE
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

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

ligand1="lig_"$(pwd | rev | cut -d "/" -f 2 | rev | cut -d "_" -f 3)
ligand2="lig_"$(pwd | rev | cut -d "/" -f 2 | rev | cut -d "_" -f 5)
i=$(pwd | rev | cut -d "/" -f 1 | rev | cut -d "_" -f 2)
j=$(pwd | rev | cut -d "/" -f 1 | rev | cut -d "_" -f 3)

# list input files
## pair pose 0 with pose 0 and pose 1 with pose 1
## but make sure that both ligands dominate the initial coordinates once for each pose pair
if (( j == 0 ))
then
file1=../../${ligand1}/ligand_premature.itp
file2=../../${ligand2}/ligand_premature.itp
file3=../../${ligand1}/pose_${i}/ligand.mol2
file4=../../${ligand2}/pose_${i}/ligand.mol2
file5=../../${ligand1}/pose_${i}/ligand.gro
file6=../../${ligand2}/pose_${i}/ligand.gro
elif (( j == 1 ))
then
file1=../../${ligand2}/ligand_premature.itp
file2=../../${ligand1}/ligand_premature.itp
file3=../../${ligand2}/pose_${i}/ligand.mol2
file4=../../${ligand1}/pose_${i}/ligand.mol2
file5=../../${ligand2}/pose_${i}/ligand.gro
file6=../../${ligand1}/pose_${i}/ligand.gro
fi

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
# For some examples, harmless atom name mismatches were observed, too
# RDKit ignores H atoms such that the alignment of ligand B on ligand A can lead to long bonds between the heavy atoms of ligand B and its H atoms
gmx grompp -f $PATH_TO_SCRIPTS/mdp/fixStructureOfHybridLigand.mdp -c merged_old.gro -r merged_old.gro -p topol_ligandInWater.top -o merged.tpr -maxwarn 3
gmx mdrun -deffnm merged -ntmpi 1 || true

if ! [ -f merged.gro ]
then
echo "For pose_${i}_${j} of edge_${ligand1}_${ligand2}, the very short energy minimisation of the newly formed hybrid ligand failed. Deleting the folder as this hybrid ligand pose should not be handled by the current workflow."
cd ..
rm -rf pose_${i}_${j}
continue
fi

# make sure that the ligand doesn't jump across periodic boundary conditions during the short EM
gmx trjconv -f merged.gro -s merged.tpr -o nojump.gro -pbc nojump <<-eof
System
eof
mv nojump.gro merged.gro

# clean up
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

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload openbabel/3.1.1
module unload ambertools/22
module unload acpype/2022.7.21
module unload gromacs/2023.2
module unload rdkit/2023.3.1

export GMXLIB=""
echo "CADD: Create Hybrid Ligands 2 -- exiting ($(pwd))" >>$LOGFILE
