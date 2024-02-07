#!/bin/bash
set -e

PATH_TO_SCRIPTS=CADD/scripts

# start simulation
if [ -f readyForAnalysisLigand]
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_ligand.sh
slurmJobId=$(ls | grep "file_ligand" | cut -d "_" -f 3)
else
sbatch ${PATH_TO_SCRIPTS}/runAWH_ligand.sh
slurmJobId=$(ls | grep "file_ligand" | cut -d "_" -f 3)
fi

checkPeriodLigand=5000 # calculate new free-energy estimate every 5 ns; TODO: turn into user parameter

# wait 5 min to be sure to stop the simulation after a checkpoint file has been written
sleep 300

while true; do

# check simulation progress once every 30 minutes
sleep 1800
numberOfValuesExpectedLigand=0

# start analysing the free-energy estimate once the respective replica leaves the initial stage of the AWH simulation
if grep -q "out of the initial stage" production_ligand.log
then
ligandOutOfInitialStage=$(grep "out of the initial stage" production_ligand.log | cut -d "=" -f 2 | cut -d "." -f 1)

lineTimeLigand=$(($(grep -n "Step           Time" production_ligand.log | tail -1 | cut -d ":" -f 1)+1 ))
timeLigand=$(head -${lineTimeLigand} production_ligand.log | tail -1 | awk '{print $2}')

# Calculate time difference between the current time and the time when AWH left the initial stage and divide by the frequency at which free-energy estimates should be calculated (interested in ceil result)
numberOfValuesExpectedLigand=$(echo "${timeLigand}" "${ligandOutOfInitialStage}" | awk '{print $1-$2}')
numberOfValuesExpectedLigand=$(echo "${numberOfValuesExpectedLigand}" "${checkPeriodLigand}" | awk '{print $1/$2}')
numberOfValuesExpectedLigand=$(( ${numberOfValuesExpectedLigand%.*}+1 ))
fi

# Stop the simulation when we should analyse
if (( $(cat ../freeEnergySummaryLigand.txt | wc -l)<${numberOfValuesExpectedLigand} ))
then
touch readyForAnalysisLigand
scancel ${slurmJobId}
break
fi

done
