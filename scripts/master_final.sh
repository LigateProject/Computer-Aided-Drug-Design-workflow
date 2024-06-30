#!/bin/bash
set -e; shopt -s expand_aliases

source ./leonardo.sh

echo "---------------------------" > $LOGFILE
echo " CADD Workflow: target=$CADD_TARGET " >> $LOGFILE
date >> $LOGFILE
hostname >> $LOGFILE
echo "---------------------------" >>$LOGFILE

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}
export CADD_SUBMISSION_DIR=$(pwd)

target=${CADD_TARGET}

# This script is for illustration purposes and was not tested yet!

step1=$(sbatch --parsable -J "1.CheckProtein" ${PATH_TO_SCRIPTS}/launchCheckProtein.sh )

###################################
# Docking with LiGen (not public) #
###################################

# Executes GROMACS_LiGen_integration_final.sh
# All these scripts can be run on 1 CPU core.
step2=$(sbatch -N1 -n1 -t 1:00:00 -p dcgp_usr_prod -A cin_staff -J "2. RIFG" --dependency=afterok:"$step1" --parsable ${PATH_TO_SCRIPTS}/runInputFileGeneration_final.sh)

# Pair ligands
step3=$(sbatch -N1 -n1 --parsable -p dcgp_usr_prod -A cin_staff -t 4:00:00 -N1 -n1 --dependency=afterok:"$step2" -J "3. HybridLigs" --wrap="${PATH_TO_SCRIPTS}/createHybridLigands.sh" )


# Executes solvate.sh, addIons.sh
# All these scripts can be run on 1 CPU core.
step4=$(sbatch -N1 -n1 -A cin_staff --parsable -p dcgp_usr_prod --dependency=afterany:"$step3" -J "4. RIFG1" "${PATH_TO_SCRIPTS}/runInputFileGeneration1.sh"  )


# Perform energy minimisation
if ! [ -d ${OUTPUT_PATH}/${target} ]
then
exit 0
fi
step5=$(sbatch --parsable -N1 -n1 -A cin_staff -t 1:00:00 -p dcgp_usr_prod --dependency=afterok:"$step4" -J "5. Minimise" --wrap="bash ${PATH_TO_SCRIPTS}/energyMinimisation.sh")

# Executes prepareEquilibration.sh
if ! [ -d ${OUTPUT_PATH}/${target} ]
then
exit 0
fi
step6=$(sbatch -N1 -n1 -A cin_staff --parsable --dependency=afterok:"$step5" -J "6. RIFG2" ${PATH_TO_SCRIPTS}/runInputFileGeneration2.sh )

# Perform equilibration
# The script can be run directly with the output of runInputFileGeneration.sh
# In principle, starting it after the SLURM job for the previous script returned with status OK should be fine. The only caveat is that error catching was implemented in the previous script such that the script does not fail when encountering expected workflow errors. Instead, it removes the directory of the affected replica and/or ligand pair and returns with exit status 0. To avoid running on an empty target directory, a check for the existence of folders named "edge_*" has to be executed.
# equilibrate.sh is expected to be run on the head node, it submits all equilibration runs needed for all replicas of all ligand pairs for the protein target under study
#if ! [ -d ${OUTPUT_PATH}/${target}/edge_* ]
#then
#exit 0
#fi
step7=$(sbatch --parsable -N1 -n1 -A cin_staff -p dcgp_usr_prod -t 8:00:00 --dependency=afterok:"$step6" -J "7.Equilibrate" --wrap="bash ${PATH_TO_SCRIPTS}/equilibrate.sh")

# Prepare the input TPR files for the final AWH simulations for all replicas of all ligand pairs of the protein target under study
# The script can be run directly with the output of equilibrate.sh
# In principle, starting it after the previous script returned with status OK should be fine. The only caveat is that error catching was implemented in the previous script such that the script does not fail when encountering expected workflow errors. Instead, it removes the directory of the affected replica and/or ligand pair and returns with exit status 0. To avoid running on an empty target directory, a check for the existence of folders named "edge_*" has to be executed.
# Requires 1 CPU core and a few minutes of execution time at most
#if ! [ -d ${OUTPUT_PATH}/${target}/edge_* ]
#then
#exit 0
#fi
step8=$(sbatch -N1 -n1 -A cin_staff --parsable  --dependency=afterok:"$step7" -J "8. PPS" ${PATH_TO_SCRIPTS}/prepareProductionSimulations.sh )

# Run AWH simulations yielding the RBFE estimate
# The script can be run directly with the output of prepareProductionSimulations.sh
# In principle, starting it after the SLURM job for the previous script returned with status OK should be fine. The only caveat is that error catching was implemented in the previous script such that the script does not fail when encountering expected workflow errors. Instead, it removes the directory of the affected replica and/or ligand pair and returns with exit status 0. To avoid running on an empty target directory, a check for the existence of folders named "edge_*" has to be executed.
# launchAWH.sh is expected to be run in the background on the head node, it submits and manages all AWH simulations needed for all replicas of all ligand pairs for the protein target under study
# The script has so far been tested for one ligand pair
#if ! [ -d ${OUTPUT_PATH}/${target}/edge_* ]
#then
#exit 0
#fi

# The sbatch command wont finish until job finished
sbatch --wait -J "9. launch AWH" -N1 -n1 -A cin_staff -p dcgp_usr_prod -t 24:00:00 --dependency=afterok:"$step8" --wrap="bash ${PATH_TO_SCRIPTS}/launchAWH.sh" &

wait  # wait for last sbatch

# Calculate ligand ranking based on RBFE estimates for all ligand pairs provided in the respective "edge_*" folder in a file called RBFE.dat
cd ${OUTPUT_PATH}/${target}
cp ${PATH_TO_SCRIPTS}/summarizeLigandRelativeFE.py .
cp ${PATH_TO_SCRIPTS}/internal_ufloat.py .
python3 summarizeLigandRelativeFE.py >> ligandRankingStdOut.txt
rm -r summarizeLigandRelativeFE.py internal_ufloat.py __pycache__

date >> $LOGFILE
echo "------ CADD WORKFLOW FINISHED -----------" >> $LOGFILE
