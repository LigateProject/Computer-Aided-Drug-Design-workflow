#!/bin/bash

#SBATCH --job-name=equilibrateLigand
#SBATCH -p boost_usr_prod
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -t 30:00
#SBATCH -A cin_staff
#SBATCH -o slurm-%J.out

set -e

pose=$(basename $PWD)
echo "CADD: Equilibrating ligand $pose" >> $LOGFILE
module use $CADD_SOFTWARE_MODULES
module load gromacs/2023.2

gmx mdrun -deffnm equiNVT_ligand
echo "CADD: Equilibrating ligand -- finished " >> $LOGFILE
