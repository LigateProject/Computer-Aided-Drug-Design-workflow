#!/bin/bash

#SBATCH --job-name=equilibrateComplex
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -t 5:00:00
#SBATCH -o slurm-%J.out

set -e

module load gromacs/2023.2 gromacs=gmx

gmx mdrun -deffnm equiNVT_complex
