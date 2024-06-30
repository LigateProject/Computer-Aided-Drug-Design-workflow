#!/bin/bash

#SBATCH --job-name=EnergyMinimisationComplex
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1 -n4
#SBATCH -G 0
#SBATCH -t 1:00:00
#SBATCH -A cin_staff
#SBATCH -o slurm-%J.out

set -e

module use $CADD_SOFTWARE_MODULES
module load gromacs/2023.2

PATH_TO_MDP=${CADD_SCRIPTS_DIR}/mdp

target=$(pwd | rev | cut -d "/" -f 1 | rev)

# run grompp to get input .tpr file
gmx grompp -f ${PATH_TO_MDP}/em_l0.mdp -c ions_complex.gro -r ions_complex.gro -p topol_amber.top -o EM_complex.tpr -po EMout_complex.mdp -maxwarn 2 || true

## error catching
## EM.tpr must exist
if ! [ -f EM_complex.tpr ]
then
echo "$(pwd): the energy minimisation failed for the complex (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
module unload gromacs/2023.2
exit 0
fi

# run energy minimisation
gmx mdrun -v -deffnm EM_complex || true

## error catching
## EM.gro must exist and must not be empty
if ! [ -s EM_complex.gro ]
then
echo "$(pwd): the energy minimisation failed for the complex (mdrun error). Deleting folder!"
cd ..
rm -rf ${pose}
module unload gromacs/2023.2
exit 0
fi

# file clean up such that subsequent scripts can be used without any problems
mv EM_complex.gro ions_complex.gro
rm EM_complex*
rm EMout_complex*

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
