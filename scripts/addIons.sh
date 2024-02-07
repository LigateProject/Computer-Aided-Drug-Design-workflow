#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=CADD/scripts/mdp

target=$(pwd | rev | cut -d "/" -f 1 | rev)

for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}
# need to add ions to get a neutral simulation box
## It should be -maxwarn 1, but right now, we have some issues with STaGE not re-naming atoms properly (ugly fix that must be replaced)
gmx grompp -f ${PATH_TO_MDP}/em_l0.mdp -c solvated_complex.gro -r solvated_complex.gro -p topol_amber.top -o addIons_complex.tpr -po mdout_complex.mdp -maxwarn 2 || true
gmx grompp -f ${PATH_TO_MDP}/em_l0.mdp -c solvated_ligand.gro -r solvated_ligand.gro -p topol_ligandInWater.top -o addIons_ligand.tpr -po mdout_ligand.mdp -maxwarn 2 || true

## error catching
## addIons.tpr must exist
if ! [ -f addIons_complex.tpr ]
then
echo "For ${pose} of ${edge}, ions could not be added to the complex structure (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -f addIons_ligand.tpr ]
then
echo "For ${pose} of ${edge}, ions could not be added to the ligand (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

set +e
gmx genion -s addIons_complex.tpr -o ions_complex.gro -p topol_amber.top -pname NA -nname CL -neutral <<-eof
SOL
eof
set -e
set +e
gmx genion -s addIons_ligand.tpr -o ions_ligand.gro -p topol_ligandInWater.top -pname NA -nname CL -neutral <<-eof
SOL
eof
set -e

## error catching
## ions.gro must exist and must not be empty
if ! [ -s ions_complex.gro ]
then
echo "For ${pose} of ${edge}, ions could not be added to the complex structure (genion error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -s ions_ligand.gro ]
then
echo "For ${pose} of ${edge}, ions could not be added to the ligand (genion error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

rm mdout_complex.mdp addIons_complex.tpr solvated_complex.gro 
rm mdout_ligand.mdp addIons_ligand.tpr solvated_ligand.gro 
if [ -f "#topol_amber.top.1#" ]
then
rm "#topol_amber.top.1#"
echo "For ${pose} of ${edge}, ions were added to the complex. Removing topology back-up."
else
echo "For ${pose} of ${edge}, no ions were added to the complex."
fi
if [ -f "#topol_ligandInWater.top.1#" ]
then
rm "#topol_ligandInWater.top.1#"
echo "For ${pose} of ${edge}, ions were added to the ligand. Removing topology back-up."
else
echo "For ${pose} of ${edge}, no ions were added to the ligand."
fi

echo "For ${pose} of ${edge}, the workflow step 'adding ions' has been completed successfully!"

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the workflow step 'adding ions' could not be completed successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "The workflow step 'adding ions' could not be completed successfully for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
