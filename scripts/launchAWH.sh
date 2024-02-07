#!/bin/bash
set -e

OUTPUT_PATH=${ligate}/CADDValidation
PATH_TO_SCRIPTS=CADD/scripts

targets=(bace_p2)

# This script starts all relevant background processes that will submit the individual AWH simulations and monitor them

cd ${OUTPUT_PATH}
for target in ${targets[@]}; do
cd ${target}
for edge in edge_*; do
cd ${edge}
touch freeEnergySummaryComplex.txt
touch freeEnergySummaryLigand.txt
for pose in pose_*; do
cd ${pose}

bash ${PATH_TO_SCRIPTS}/monitorAWH_complex.sh
echo $! >> pid_complex.txt
bash ${PATH_TO_SCRIPTS}/monitorAWH_ligand.sh
echo $! >> pid_ligand.txt

cd ..
done # poses
cd ..
done # edges
cd ..
done # targets

while true; do
# Every minute,
sleep 60
# go through all folders and check whether a free-energy can be computed
for target in ${targets[@]}; do
cd ${target}
countEdge=0
for edge in edge_*; do
cd ${edge}

if [ -f isConvergedComplex ] && [ -f isConvergedLigand ]
then
countEdge=$(( ${countEdge}+1 ))
continue
fi

countComplex=0
countLigand=0
for pose in pose_*; do
cd ${pose}

if [ -f readyForAnalysisComplex ]
then
countComplex=$(( ${countComplex}+1 ))
fi
if [ -f readyForAnalysisLigand ]
then
countLigand=$(( ${countLigand}+1 ))
fi

cd ..
done # poses

if (( ${countComplex} == $(ls | grep -c "pose") ))
then
bash ${PATH_TO_SCRIPTS}/analyseAWH_complex.sh
fi
if (( ${countLigand} == $(ls | grep -c "pose") ))
then
bash ${PATH_TO_SCRIPTS}/analyseAWH_ligand.sh
fi

if [ -f isConvergedComplex ] && [ -f isConvergedLigand ]
then
RBFE=$(echo $(tail -1 freeEnergySummaryComplex.txt | awk '{print $2}') $(tail -1 freeEnergySummaryligand.txt | awk '{print $2}') | awk '{print $1-$2}')
stderrComplex=$(echo $(tail -1 freeEnergySummaryComplex.txt) | awk '{print $3}')
stderrLigand=$(echo $(tail -1 freeEnergySummaryLigand.txt) | awk '{print $3}')
stderr=$(echo $(echo ${stderrComplex} ${stderrComplex} | awk '{print $1*$2}') $(echo ${stderrLigand} ${stderrLigand} | awk '{print $1*$2}') | awk '{print $1+$2}')
stderr=$(echo "sqrt(${s})" | bc)
echo "${RBFE} +- ${stderr}" >> RBFE.dat
fi

cd ..
done # edges

if (( ${countEdge} == $(ls | grep -c "edge") ))
then
continue
fi

cd ..
done # targets

done
