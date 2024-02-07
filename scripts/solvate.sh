#!/bin/bash
set -e

module load gromacs/2023.2 gromacs=gmx

target=$(pwd | rev | cut -d "/" -f 1 | rev)

for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

# place the complex in a rhombic dodecahedron of the desired size (1.5 nm distance to the edges)
gmx editconf -f full.gro -o correctBox_complex.gro -bt dodecahedron -d 1.5 || true
gmx editconf -f merged.gro -o correctBox_ligand.gro -bt dodecahedron -d 1.5 || true
rm full.gro merged.gro

## error catching
## correctBox.gro must exist and must not be empty
if ! [ -s correctBox_complex.gro ]
then
echo "For ${pose} of ${edge}, the complex structure could not be placed in a box of the correct size. Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -s correctBox_ligand.gro ]
then
echo "For ${pose} of ${edge}, the ligand could not be placed in a box of the correct size. Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

# solvate with TIP3P water
gmx solvate -cp correctBox_complex.gro -cs spc216.gro -p topol_amber.top -o solvated_complex.gro || true
gmx solvate -cp correctBox_ligand.gro -cs spc216.gro -p topol_ligandInWater.top -o solvated_ligand.gro || true

## error catching
## the solvated .gro file must exist and must not be empty
if ! [ -s solvated_complex.gro ]
then
echo "For ${pose}, the complex structure could not be solvated (.gro file not obtained). Deleting folder!"
cd ..
rm -rf ${pose}
continue
## water molecules must be added to the topology
elif ! $(diff 'topol_amber.top' '#topol_amber.top.1#' | grep -q "SOL")
then
echo "For ${pose} of ${edge}, the complex structure could not be solvated (water molecules were not added to topology). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
## the solvated .gro file must exist and must not be empty
if ! [ -s solvated_ligand.gro ]
then
echo "For ${pose}, the complex structure could not be solvated (.gro file not obtained). Deleting folder!"
cd ..
rm -rf ${pose}
continue
## water molecules must be added to the topology
elif ! $(diff 'topol_ligandInWater.top' '#topol_ligandInWater.top.1#' | grep -q "SOL")
then
echo "For ${pose} of ${edge}, the complex structure could not be solvated (water molecules were not added to topology). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

rm "#topol_amber.top.1#" correctBox_complex.gro
rm "#topol_ligandInWater.top.1#" correctBox_ligand.gro
sed -i solvated_complex.gro -e "s/HOH/SOL/g"
sed -i solvated_ligand.gro -e "s/HOH/SOL/g"

echo "For ${pose} of ${edge}, the complex and the ligand were solvated successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, no pose could be solvated successfully. Deleting directory!"
cd ..
rm -rf ${edge}
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "No edge could be solvated successfully. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
