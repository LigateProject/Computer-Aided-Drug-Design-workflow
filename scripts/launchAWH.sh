#!/bin/bash
set -e

echo "CADD: launch AWH " >> $LOGFILE
OUTPUT_PATH=${CADD_OUTPUT_DIR}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

target=${CADD_TARGET}

# This script starts all relevant background processes that will submit the individual AWH simulations and monitor them

cd ${OUTPUT_PATH}/${target}
for edge in edge_*; do
cd ${edge}
touch freeEnergySummaryComplex.txt
touch freeEnergySummaryLigand.txt
for pose in pose_*; do
cd ${pose}

bash ${PATH_TO_SCRIPTS}/monitorAWH_complex.sh &
echo $! >> pid_complex.txt
bash ${PATH_TO_SCRIPTS}/monitorAWH_ligand.sh &
echo $! >> pid_ligand.txt

cd ..
done # poses
cd ..
done # edges

while true; do
# Every minute,
sleep 60
# go through all folders and check whether a free-energy can be computed
countEdge=0
for edge in edge_*; do
cd ${edge}

if [ -f isConvergedComplex ] && [ -f isConvergedLigand ]
then
cd ..
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

if (( ${countComplex} == $(ls | grep -c "pose") )) && ! [ -f isConvergedComplex ]
then
echo "Calculating new free energy value for the complex ..."
bash ${PATH_TO_SCRIPTS}/analyseAWH_complex.sh
# re-start simulation if not converged
## pid_complex.txt is only removed if the simulation has converged
testPose=$(ls | grep "pose_" | head -1)
if [ -f ${testPose}/pid_complex.txt ]
then
for pose in pose*; do
cd ${pose}
rm pid_complex.txt
bash ${PATH_TO_SCRIPTS}/monitorAWH_complex.sh &
echo $! >> pid_complex.txt
cd ..
done
fi
fi

# final clean-up
if [ -f isConvergedComplex ]
then
for pose in pose_*; do
cd ${pose}
## we only want to clean up once
if ! [ -f readyForAnalysisComplex ]
then
cd ..
break
fi
rm readyForAnalysisComplex
rm -r cpt_complex
cd ..
done
fi

if (( ${countLigand} == $(ls | grep -c "pose") )) && ! [ -f isConvergedLigand ]
then
echo "Calculating new free energy value for the ligand ..."
bash ${PATH_TO_SCRIPTS}/analyseAWH_ligand.sh
# re-start simulation if not converged
## pid_ligand.txt is only removed if the simulation has converged
testPose=$(ls | grep "pose_" | head -1)
if [ -f ${testPose}/pid_ligand.txt ]
then
for pose in pose*; do
cd ${pose}
rm pid_ligand.txt
bash ${PATH_TO_SCRIPTS}/monitorAWH_ligand.sh &
echo $! >> pid_ligand.txt
cd ..
done
fi
fi

# final clean-up
if [ -f isConvergedLigand ]
then
for pose in pose_*; do
cd ${pose}
## we only want to clean up once
if ! [ -f readyForAnalysisLigand ]
then
cd ..
break
fi
rm readyForAnalysisLigand
rm -r cpt_ligand
cd ..
done
fi

if [ -f isConvergedComplex ] && [ -f isConvergedLigand ]
then
RBFE=$(echo $(tail -1 freeEnergySummaryComplex.txt | awk '{print $2}') $(tail -1 freeEnergySummaryLigand.txt | awk '{print $2}') | awk '{print $1-$2}')
stderrComplex=$(echo $(tail -1 freeEnergySummaryComplex.txt) | awk '{print $3}')
stderrLigand=$(echo $(tail -1 freeEnergySummaryLigand.txt) | awk '{print $3}')
stderr=$(echo $(echo ${stderrComplex} ${stderrComplex} | awk '{print $1*$2}') $(echo ${stderrLigand} ${stderrLigand} | awk '{print $1*$2}') | awk '{print $1+$2}')
stderr=$(echo "sqrt(${stderr})" | bc)
echo "${RBFE} +- ${stderr}" >> RBFE.dat
fi

cd ..
done # edges

# stop script if all RBFE calculations have converged
if (( ${countEdge}==$(ls | grep -c "edge")  ))
then
echo "RBFE calculations for all edges have converged!"
exit
fi

done
echo "CADD: launch AWH --exiting " >> $LOGFILE
