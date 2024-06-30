#!/bin/bash
set -e

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

# remove remainders from a previous script execution
if (( $(ls file_ligand_* | wc -w) > 0 ))
then
rm file_ligand_*
fi

# start simulation
if [ -f readyForAnalysisLigand ]
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_ligand.sh
rm readyForAnalysisLigand
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_ligand_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_ligand" | cut -d "_" -f 3)
# remove old CPT files but keep the ten most recent for re-starts
if (( $(ls -lt cpt_ligand/ | wc -l) > 11 ))
then
cd cpt_ligand
rm $(ls -lt | tail -n +12 | rev | cut -d " " -f 1 | rev)
cd ..
fi
else
sbatch ${PATH_TO_SCRIPTS}/runAWH_ligand.sh
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_ligand_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_ligand" | cut -d "_" -f 3)
fi

# calculate new free-energy estimate every 2 ns; TODO: turn into user parameter
# must be a multiple of awh-nstout (currently 100 ps)
checkPeriodLigand=2000

while true; do

# check simulation progress every minute
sleep 60
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
echo "Interrupting ligand replica for free energy calculation (${slurmJobIdToBeCancelled}, $(pwd))!"
touch readyForAnalysisLigand
scancel ${slurmJobIdToBeCancelled}
break
# Re-start if simulation is no longer running or in the queue but should be
else
if ! $(squeue | grep -q "${slurmJobIdToBeCancelled}")
then
echo "Re-starting failed ligand simulation (${slurmJobIdToBeCancelled}, $(pwd))!"
# remove remainders from previous run
if (( $(ls file_ligand_* | wc -w) > 0 ))
then
#cat file_ligand_*
rm file_ligand_*
fi
# if we re-start, let's be cautious and assume that the last CPT file is corrupted
# in worst case, we delete the state of a long simulation and have to re-start from scratch
if $(ls -lt cpt_complex | grep -q "state")
then
cd cpt_ligand
rm $(ls -ltr | tail -1 | rev | cut -d " " -f 1 | rev)
cd ..
fi
# remove old CPT files but keep the ten most recent for re-starts
if (( $(ls -lt cpt_ligand/ | wc -l) > 11 ))
then
cd cpt_ligand
rm $(ls -lt | tail -n +12 | rev | cut -d " " -f 1 | rev)
cd ..
fi
# re-start
if $(ls -lt cpt_ligand | grep -q "state")
then
sbatch ${PATH_TO_SCRIPTS}/restartAWH_ligand.sh
else
for targetFile in "production_ligand.edr" "production_ligand.log" "production_ligand.xtc" "production_ligand.xvg"; do
if [ -f ${targetFile} ]
then
rm ${targetFile}
fi
done
sbatch ${PATH_TO_SCRIPTS}/runAWH_ligand.sh
fi
# wait until SLURM job has started
while true; do
sleep 10
if [ -f file_ligand_* ]
then
break
fi
done
slurmJobIdToBeCancelled=$(ls | grep "file_ligand" | cut -d "_" -f 3)
fi
fi

done
