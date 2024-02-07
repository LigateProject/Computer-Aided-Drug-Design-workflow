#!/bin/bash
set -e

PATH_TO_SCRIPTS=CADD/scripts

# start simulation
if [ -f readyForAnalysisComplex]
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_complex.sh
slurmJobId=$(ls | grep "file_complex" | cut -d "_" -f 3)
else
sbatch ${PATH_TO_SCRIPTS}/runAWH_complex.sh
slurmJobId=$(ls | grep "file_complex" | cut -d "_" -f 3)
fi

checkPeriodComplex=10000 # calculate new free-energy estimate every 10 ns; TODO: turn into user parameter

# wait 5 min to be sure to stop the simulation after a checkpoint file has been written
sleep 300

while true; do

# check simulation progress once every 60 minutes
sleep 3000
numberOfValuesExpectedComplex=0

# start analysing the free-energy estimate once the respective replica leaves the initial stage of the AWH simulation
if grep -q "out of the initial stage" production_complex.log
then
complexOutOfInitialStage=$(grep "out of the initial stage" production_complex.log | cut -d "=" -f 2 | cut -d "." -f 1)

lineTimeComplex=$(($(grep -n "Step           Time" production_complex.log | tail -1 | cut -d ":" -f 1)+1 ))
timeComplex=$(head -${lineTimeComplex} production_complex.log | tail -1 | awk '{print $2}')

# Calculate time difference between the current time and the time when AWH left the initial stage and divide by the frequency at which free-energy estimates should be calculated (interested in ceil result)
numberOfValuesExpectedComplex=$(echo "${timeComplex}" "${complexOutOfInitialStage}" | awk '{print $1-$2}')
numberOfValuesExpectedComplex=$(echo "${numberOfValuesExpectedComplex}" "${checkPeriodComplex}" | awk '{print $1/$2}')
numberOfValuesExpectedComplex=$(( ${numberOfValuesExpectedComplex%.*}+1 ))
fi

# Stop the simulation when we should analyse
if (( $(cat ../freeEnergySummaryComplex.txt | wc -l)<${numberOfValuesExpectedComplex} ))
then
touch readyForAnalysisComplex
scancel ${slurmJobId}
break
fi

done
