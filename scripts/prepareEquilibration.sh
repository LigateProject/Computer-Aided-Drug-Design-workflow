#!/bin/bash
set -e

export LC_NUMERIC="en_US.UTF-8"

module load gromacs/2023.2 gromacs=gmx

PATH_TO_MDP=CADD/scripts/mdp

target=$(pwd | rev | cut -d "/" -f 1 | rev)

for edge in edge_*; do
cd ${edge}
for pose in pose_*; do
cd ${pose}

# check whether the alchemical transformation alters the charge (cannot be more than one unit) and add a buffer particle if so
## difference chargeA - chargeB
number1=$(( $(grep -n "\[ atoms \]" merged.itp | cut -d ":" -f 1)+2 ))
number2=$(( $(grep -n "\[ bonds \]" merged.itp | cut -d ":" -f 1)-2 ))

chargeA=($(head -${number2} merged.itp | tail -n +${number1} | awk '{print $7}'))
chargeB=($(head -${number2} merged.itp | tail -n +${number1} | awk '{print $10}'))

sumChargeA=0
sumChargeB=0

for index in $(seq 0 $(( ${#chargeA[@]}-1 ))); do
sumChargeA=$(echo ${sumChargeA} ${chargeA[${index}]} | awk '{print $1+$2}')
sumChargeB=$(echo ${sumChargeB} ${chargeB[${index}]} | awk '{print $1+$2}')
done

sumChargeA=$(printf "%.0f" "${sumChargeA}")
sumChargeB=$(printf "%.0f" "${sumChargeB}")

chargeDifference=$(echo ${sumChargeA} ${sumChargeB} | awk '{$1-$2}')
chargeDifference=$(printf "%.0f" "${chargeDifference}")

if (( ${chargeDifference} != 0 ))
then
## create buffer particle to compensate charge difference
### add buffer particle to ffMOL.itp
echo " BUF      BUF         0.00000  0.00000   A     0.00000e+00   0.00000e+00" >> ffMOL.itp

### add buffer particle to ligand topology (no bonded interactions, no mass, no vdW, just charge)
lineNumber=$(( $(grep -n "\[ bonds \]" merged.itp | cut -d ":" -f 1)-2 ))
lastLine=$(head -${lineNumber} merged.itp | tail -1)
element1=$(( $(echo ${lastLine} | awk '{print $1}')+1 ))
element2="BUF"
element3=$(echo ${lastLine} | awk '{print $3}')
element4=$(echo ${lastLine} | awk '{print $4}')
element5="BUF"
element6=$(( $(echo ${lastLine} | awk '{print $6}')+1 ))
element7=0 # charge
element8=1 # mass
element9="BUF"
element10=${chargeDifference} # charge
element11=1 # mass

newLine=$(printf "%6d%12s%7d%7s%7s%7d%11.6f%11.5f%12s%11.6f%11.5f\n" ${element1} ${element2} ${element3} ${element4} ${element5} ${element6} ${element7} ${element8} ${element9} ${element10} ${element11})
head -${lineNumber} merged.itp >> intermediary.itp
echo "${newLine}" >> intermediary.itp
tail -n +$(( ${lineNumber}+1 )) merged.itp >> intermediary.itp
mv intermediary.itp merged.itp

### add buffer particle to .gro files
#### approach: add half the box dimension to the position of the last ligand atom
for gro in ions_complex.gro ions_ligand.gro; do
x=$(tail -1 ${gro} | awk '{print $1}')
y=$(tail -1 ${gro} | awk '{print $2}')
z=$(tail -1 ${gro} | awk '{print $3}')

x_half=$(echo ${x} 2 | awk '{print $1/$2}')
y_half=$(echo ${y} 2 | awk '{print $1/$2}')
z_half=$(echo ${z} 2 | awk '{print $1/$2}')

lineNumber=$(grep -n "LIG" ${gro} | tail -1 | cut -d ":" -f 1)
oldLine=$(head -${lineNumber} ${gro} | tail -1)
element1=$(echo ${oldLine} | awk '{print $1}')
element2="BUF"
element3=$(( $(echo ${oldLine} | awk '{print $3}')+1 ))
element4=$(echo ${oldLine} | awk '{print $4}')
element4=$(echo ${element4} ${x_half} | awk '{print $1+$2}')
element5=$(echo ${oldLine} | awk '{print $5}')
element5=$(echo ${element5} ${y_half} | awk '{print $1+$2}')
element6=$(echo ${oldLine} | awk '{print $6}')
element6=$(echo ${element6} ${z_half} | awk '{print $1+$2}')
newLine=$(printf "%8s%7s%5d%8.3f%8.3f%8.3f\n" ${element1} ${element2} ${element3} ${element4} ${element5} ${element6})

head -${lineNumber} ${gro} >> intermediary.gro
echo "${newLine}" >> intermediary.gro
tail -n +$(( ${lineNumber}+1 )) ${gro} >> intermediary.gro
mv intermediary.gro ${gro}

#### update the number of atoms
head -1 ${gro} >> intermediary.gro
newLine=$(printf "%5d\n" $(( $(head -2 ${gro} | tail -1 | awk '{print $1}')+1 )))
echo "${newLine}" >> intermediary.gro
tail -n +3 ${gro} >> intermediary.gro
mv intermediary.gro ${gro}

#### let GROMACS renumber
gmx editconf -f ${gro} -o intermediary.gro
mv intermediary.gro ${gro}
done
fi

# run grompp to get input .tpr file
# need to use soft-core potentials although vdW is not switched off => ignore the corresponding warning
# system may have a slight net charge due to rounding errors => ignore that warning, too
gmx grompp -f ${PATH_TO_MDP}/eq_nvt_l0.mdp -c ions_complex.gro -p topol_amber.top -o equiNVT_complex.tpr -po equiNVTOut_complex.mdp -maxwarn 2 || true
gmx grompp -f ${PATH_TO_MDP}/eq_nvt_l0.mdp -c ions_ligand.gro -p topol_ligandInWater.top -o equiNVT_ligand.tpr -po equiNVTOut_ligand.mdp -maxwarn 2 || true
## error catching
## TPR file for equilibration must exist
if ! [ -f equiNVT_complex.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the complex (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi
if ! [ -f equiNVT_ligand.tpr ]
then
echo "For ${pose} of ${edge}, the TPR file could not be generated for the ligand (grompp error). Deleting folder!"
cd ..
rm -rf ${pose}
continue
fi

rm equiNVTOut_complex.mdp
rm equiNVTOut_ligand.mdp

echo "For ${pose} of ${edge}, the TPR file has been generated successfully."

cd ..
done # poses

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) == 0 ));
then
echo "For ${edge}, the TPR file could not be generated successfully for any of the poses. Deleting directory!"
cd ..
rm -rf ${edge}
fi

cd ..
done # edges

# remove directory if none of the edges survives
if (( $(ls -lh | grep "edge_" | wc -l) == 0 ));
then
echo "The TPR file could not be generated successfully for any of the edges. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# unload required external software to restore the environment at the beginning of the script
module unload gromacs/2023.2
