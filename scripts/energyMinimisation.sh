#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=CADD/scripts/mdp

target=$(pwd | rev | cut -d "/" -f 1 | rev)

for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}
# run grompp to get input .tpr file
gmx grompp -f ${PATH_TO_MDP}/em_l0.mdp -c ions_complex.gro -r ions_complex.gro -p topol_amber.top -o EM_complex.tpr -po EMout_complex.mdp -maxwarn 2 || true
gmx grompp -f ${PATH_TO_MDP}/em_l0.mdp -c ions_ligand.gro -r ions_ligand.gro -p topol_ligandInWater.top -o EM_ligand.tpr -po EMout_ligand.mdp -maxwarn 2 || true

## error catching
## EM.tpr must exist
if ! [ -f EM_complex.tpr ]
then
echo "For ${pose} of ${edge}, the energy minimisation failed for the complex (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -f EM_ligand.tpr ]
then
echo "For ${pose} of ${edge}, the energy minimisation failed for the ligand (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

# run energy minimisation
gmx mdrun -v -deffnm EM_complex || true
gmx mdrun -v -deffnm EM_ligand || true

## error catching
## EM.gro must exist and must not be empty
if ! [ -s EM_complex.gro ]
then
echo "For ${pose} of ${edge}, the energy minimisation failed for the complex (mdrun error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -s EM_ligand.gro ]
then
echo "For ${pose} of ${edge}, the energy minimisation failed for the ligand (mdrun error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

# file clean up such that subsequent scripts can be used without any problems
mv EM_complex.gro ions_complex.gro
mv EM_ligand.gro ions_ligand.gro
rm EM*

echo "For ${edge}, the energy minimisation has been completed successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the energy minimisation could not be completed successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "The energy minimisation could not be completed successfully for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
