#!/bin/bash

#SBATCH --job-name=runAWH_ligand
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -t 48:00:00
#SBATCH -o slurm-%J.out

set -e

module load gromacs/2023.2 gromacs=gmx

touch file_ligand_${SLURM_JOB_ID}
mkdir -p cpt_ligand
gmx mdrun -deffnm production_ligand -cpnum -cpt 30 -cpo cpt_ligand/state -pin on -ntomp ${SLURM_JOB_CPUS_PER_NODE} -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu
