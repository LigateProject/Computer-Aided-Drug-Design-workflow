#!/bin/bash

#SBATCH --job-name=runAWH_complex
#SBATCH -N 1
#SBATCH -G 1
#SBATCH -t 48:00:00
#SBATCH -o slurm-%J.out

set -e

module load gromacs/2023.2 gromacs=gmx

rm file_complex*
touch file_complex_${SLURM_JOB_ID}
CPT_FILE=$(ls -ltr cpt_complex | tail -1)
CPT_FILE=$(echo ${CPT_FILE} | awk '{print $NF}')
gmx mdrun -deffnm production_complex -cpnum -cpt 30 -cpo cpt_complex/state -pin on -ntomp ${SLURM_JOB_CPUS_PER_NODE} -nb gpu -pme gpu -pmefft gpu -bonded gpu -update gpu -cpi cpt_complex/${CPT_FILE}
