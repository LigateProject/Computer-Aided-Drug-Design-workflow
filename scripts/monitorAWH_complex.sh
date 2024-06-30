#!/bin/bash
set -e

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

# remove remainders from a previous script execution
if (( $(ls file_complex_* | wc -w) > 0 ))
then
rm file_complex_*
fi

# start simulation
if [ -f readyForAnalysisComplex ]
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_complex.sh
rm readyForAnalysisComplex
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_complex_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_complex" | cut -d "_" -f 3)
# remove old CPT files but keep the ten most recent for re-starts
if (( $(ls -lt cpt_complex/ | wc -l) > 11 ))
then
cd cpt_complex
rm $(ls -lt | tail -n +12 | rev | cut -d " " -f 1 | rev)
cd ..
fi
else
sbatch ${PATH_TO_SCRIPTS}/runAWH_complex.sh
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_complex_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_complex" | cut -d "_" -f 3)
fi

# calculate new free-energy estimate every 5 ns; TODO: turn into user parameter
# must be a multiple of awh-nstout (currently 100 ps)
checkPeriodComplex=5000

while true; do

# check simulation progress every minute
sleep 60
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
echo "Interrupting complex replica for free energy calculation (${slurmJobIdToBeCancelled}, $(pwd))!"
touch readyForAnalysisComplex
scancel ${slurmJobIdToBeCancelled}
break
# Re-start if simulation is no longer running or in the queue but should be
else
if ! $(squeue | grep -q "${slurmJobIdToBeCancelled}")
then
echo "Re-starting failed complex simulation (${slurmJobIdToBeCancelled}, $(pwd))!"
# remove remainders from previous run
if (( $(ls file_complex_* | wc -w) > 0 ))
then
#cat file_complex_*
rm file_complex_*
fi
# if we re-start, let's be cautious and assume that the last CPT file is corrupted
# in worst case, we delete the state of a long simulation and have to re-start from scratch
if $(ls -lt cpt_complex | grep -q "state")
then
cd cpt_complex
rm $(ls -ltr | tail -1 | rev | cut -d " " -f 1 | rev)
cd ..
fi
# remove old CPT files but keep the ten most recent for re-starts
if (( $(ls -lt cpt_complex/ | wc -l) > 11 ))
then
cd cpt_complex
rm $(ls -lt | tail -n +12 | rev | cut -d " " -f 1 | rev)
cd ..
fi
# re-start
if $(ls -lt cpt_complex | grep -q "state")
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_complex.sh
else
for targetFile in "production_complex.edr" "production_complex.log" "production_complex.xtc" "production_complex.xvg"; do
if [ -f ${targetFile} ]
then
rm ${targetFile}
fi
done
sbatch ${PATH_TO_SCRIPTS}/runAWH_complex.sh
fi
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_complex_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_complex" | cut -d "_" -f 3)
fi
fi

done
