#!/bin/bash

#SBATCH --job-name=runAWH_complex
#SBATCH -p boost_usr_prod
#SBATCH -N 1 -n1 -c8
#SBATCH -G 1
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out
#SBATCH -A lig8_md_1

set -e

echo "CADD: restart AWH complex started" >> $LOGFILE
module use $CADD_SOFTWARE_MODULES
module load gromacs/2023.2

touch file_complex_${SLURM_JOB_ID}
#hostname >> file_complex_${SLURM_JOB_ID}
CPT_FILE=$(ls -ltr cpt_complex | tail -1)
CPT_FILE=$(echo ${CPT_FILE} | awk '{print $NF}')
# A checkpoint file is needed to restart an AWH simulation that was interrupted to calculate the RBFE estimate.
# Given the frequency at which RBFE estimates are calculated and the performance to be expected for a simulation of a protein-hybrid ligand complex, a new checkpoint file should be written every 1 min.
# TODO: turn the time interval after which checkpoint files are written (currently set to 1 min: -cpt 1) into a user parameter because it is coupled to the frequency at which free energy estimates are calculated
gmx mdrun -deffnm production_complex -cpnum -cpt 1 -cpo cpt_complex/state -pin on -ntmpi 1 -ntomp ${SLURM_JOB_CPUS_PER_NODE} -nb gpu -pme gpu -pmefft gpu -bonded gpu -update cpu -cpi cpt_complex/${CPT_FILE}
echo "CADD: restart AWH complex finished" >> $LOGFILE
