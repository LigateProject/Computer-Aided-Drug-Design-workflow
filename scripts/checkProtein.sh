#!/bin/bash

#SBATCH --job-name=checkTarget
#SBATCH -p dcgp_usr_prod
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -G 0
#SBATCH --mem=8GB
#SBATCH -t 24:00:00
#SBATCH -o slurm-%J.out

set -e

echo "CADD: CheckProtein " >>$LOGFILE
# path to additional modules
module use ${CADD_SOFTWARE_MODULES}

# make sure we start from a clean environment (mainly needed if the previous complex was deleted)
modules=(biopython/1.81 boost/1.82.0 acpype/2022.7.21 ambertools/22 openbabel/3.1.1 gromacs/2023.2 openmm/7.7.0 ost/2.4 pdb-tools/2.5.0 promod/3.3 python/3.11.1 rmsd/1.5.1 stage/1.0.0 tmbed/1.0.0)
for mod in ${modules[@]}; do
if $(module list | grep -q "${mod}")
then
module unload ${mod}
fi
done

# required external software
module load python/3.11.1
module load biopython/1.81
module load tmbed/1.0.0
module load boost/1.82.0
module load openmm/7.7.0
module load ost/2.4
module load promod/3.3
module load gromacs/2023.2
module load rmsd/1.5.1
module load pdb-tools/2.5.0

module load cuda
export OPENMM_CUDA_COMPILER=$(which nvcc)

target=$(pwd | rev | cut -d "/" -f 1 | rev)

PATH_TO_SLURM=${CADD_SUBMISSION_DIR}
INPUT_PATH=${CADD_INPUT_DIR}/${CADD_TARGET}
PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

echo "Preparing protein ..." >> $LOGFILE
file=${INPUT_PATH}/protein.pdb

canonicalFastaFile=false
# For now, we're expecting users to provide a canonical fasta file as canonical.fasta
if [ -f ${INPUT_PATH}/canonical.fasta ]
then
canonicalFastaFile=true
fi

## Sometimes, entire chains aren't modelled in the PDB file.
## Let's assume that these chains aren't important for ligand binding in this case (If a chain isn't modelled at all, we have no checks because we have no idea where it is.)
## Therefore, let's generate FASTA files from the PDB files first and compare them to the official FASTA files in the PDB later.
# Create FASTA file
cp ${file} protein.pdb
python3 ${PATH_TO_SCRIPTS}/pdb2fasta.py protein.pdb

# clean up folder if non-standard residues are found
if [ -f nonStandardResidueFound ]
then
echo "Non-standard amino acids were found in the structure. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

## We only want to remove a protein from the data set if the ligand-binding core is a membrane protein
## We don't want to remove it due to membrane anchors etc. irrelevant to the protein-ligand interaction
## Therefore, we use our self-made FASTA file for this check
# Check whether it is a membrane protein
tmbed predict -f protein.fasta -p protein.pred --out-format 0 --no-use-gpu

prediction=$(tail -1 protein.pred)
if echo ${prediction} | grep -q "B" || echo ${prediction} | grep -q "b" || echo ${prediction} | grep -q "H" || echo ${prediction} | grep -q "h" || echo ${prediction} | grep -q "S"
then
echo "Some residues of the target protein should be located in a membrane. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

## We have a data set of proteins co-crystallised with ligands.
## Therefore, N- or C-terminal extensions folding back and binding to the ligand should be visible in the structure
## The most probable scenario is that these extensions are either disordered or membrane anchors/tethers/etc., which we can't handle with our current workflow anyways
## Consequently, we only need to model missing residues inside the ligand-binding core unit

if ! ${canonicalFastaFile}
then
echo "No canonical FASTA file found. Cannot repair or close gaps in protein structure. If you do not provide a canonical FASTA file, you must provide a protein structure without gaps in which all relevant atoms except hydrogens have been resolved."
fi

# Collect information about gaps in the structure
pdb_delhetatm protein.pdb > pdb_nohetatm.pdb
pdb_gap pdb_nohetatm.pdb > gaps.txt

if ! ${canonicalFastaFile}
then
if [[ $(tail -1 gaps.txt) != "Found 0 gap(s) in the structure" ]]
then
echo "Input protein structure has gaps, but no canonical FASTA file was provided. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# clean up
rm gaps.txt pdb_nohetatm.pdb protein.fasta protein.pred
touch missingResidues.txt missingAtoms.txt residuesNotRemodelled.txt CONECT.txt

elif ${canonicalFastaFile}
then
# add missing residues to sequence of core unit
python3 ${PATH_TO_SCRIPTS}/extractSequence.py protein.fasta canonical.fasta protein.pdb pdb_nohetatm.pdb gaps.txt
rm gaps.txt
rm pdb_nohetatm.pdb

# delete folder if the official PDB FASTA sequence contains the letter "X"
if [ -f xInFastaSequence ]
then
cd ..
rm -rf ${target}
echo "The canonical FASTA sequence does not specify the identity of some residues in the sequence (X as one letter code). Stopping workflow execution!"
exit 0
fi

# delete folder if the sequence of the PDB structure cannot be properly aligned to the full sequence
if [ -f alignmentError ]
then
cd ..
rm -rf ${target}
echo "The check whether residues are missing could not be executed correctly. Stopping workflow execution!"
exit 0
fi

mv new.fasta protein.fasta

## We want to stay as close to the experimental structure as possible, i.e. do not remove crystal waters and ions (ligands are in a different file)
## We still have to run the protocol all the time because many structures miss just a few atoms in some side chains
## It seems that promod3 introduces very little changes to the experimental structure with the current settings (even though it rebuilds the side chains)
## We assume that the changes are so small that even side chain atoms complexated by metal ions do not move enough for the complexation to become impossible
## Strategy: accept current promod3 script, copy promod3 result back to PDB file (it's very unlikely to have ions or crystal waters resolved next to a missing loop), check RMSD to experimental structure
# fix protein structure
echo "CADD: Starting promod3" >> $LOGFILE
python3 ${PATH_TO_SCRIPTS}/structureNormalisation.py protein.pdb protein.fasta protein_normalised.pdb

# delete folder if the sequence of the PDB structure cannot be properly aligned to the full sequence
if [ -f normalisationError ]
then
cd ..
rm -rf ${target}
echo "CADD: promod3 didn't work correctly (alignment error). Stopping workflow execution!" >> $LOGFILE
echo "promod3 didn't work correctly (alignment error). Stopping workflow execution!"
exit 0
fi

# delete folder if input PDB file can't be read, e.g., because some atoms have improper naming in the original PDB
if [ -f normalisationErrorReadingPDB ]
then
cd ..
rm -rf ${target}
echo "promod3 didn't work correctly (error reading input PDB file). Stopping workflow execution!"
exit 0
fi

# get missing residues from promod output (more reliable than tmbed)
tailNumber=$(grep -n "Starting promod3" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
ringPunchNumber=${tailNumber}
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "ring punch(es)")
then
ringPunchNumber=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -n "ring punch(es)" | cut -d ":" -f 1)
fi
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -q "Resolved")
then
gaps=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -c "Resolved")
fi
touch missingResidues.txt

for i in $(seq 1 ${gaps}); do
modelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Resolved" | tail -${i} | head -1)
chain=$(echo ${modelled} | cut -d "." -f 1)
chain=${chain: -1}
lowerBound=$(echo ${modelled} | cut -d "." -f 2 | cut -d "-" -f 1 | tr -d -c 0-9) # promod3 uses consecutive numbering without letters
upperBound=$(echo ${modelled} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9) # promod3 uses consecutive numbering without letters
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep -q "Merged gap")
then
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep -q " ${chain}")
then
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep -q ${lowerBound})
then
lowerBound="${lowerBound} "$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep ${lowerBound} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9)
upperBound=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | head -${ringPunchNumber} | grep "Merged gap" | grep " ${chain}" | grep ${upperBound} | cut -d "." -f 3 | cut -d " " -f 1 | tr -d -c 0-9)" ${upperBound}"
fi
fi
fi
for h in $(seq 1 $(echo ${lowerBound} | wc -w)); do
for j in $(seq $(( $(echo ${lowerBound} | cut -d " " -f ${h}) + 1 )) $(( $(echo ${upperBound} | cut -d " " -f ${h}) - 1 ))); do
if (( ${j} == ($(echo ${upperBound} | cut -d " " -f ${h}) - 1) )) && (( ${i} == ${gaps} )) && (( ${h} == $(echo ${lowerBound} | wc -w) ))
then
echo -n "(group Protein and chain ${chain} and resnr ${j})" >> missingResidues.txt
else
echo -n "(group Protein and chain ${chain} and resnr ${j}) or " >> missingResidues.txt
fi
done
done
done

# get residues that promod3 refuses to remodel
touch residuesNotRemodelled.txt
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "incomplete backbone")
then
numberNotRemodelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "incomplete backbone")
for i in $(seq 1 ${numberNotRemodelled}); do
notRemodelled=$(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep "incomplete backbone" | head -${i} | tail -1)
currentResidueNumber=$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1)
currentResidueNumber=${currentResidueNumber:3}
if (( ${i} == ${numberNotRemodelled} ))
then
echo -n "(group Protein and chain "$(echo ${notRemodelled} | cut -d "." -f 1)" and resname "$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1 | cut -c1-3)" and resnr "${currentResidueNumber}")" >> residuesNotRemodelled.txt
else
echo -n "(group Protein and chain "$(echo ${notRemodelled} | cut -d "." -f 1)" and resname "$(echo ${notRemodelled} | cut -d "." -f 2 | cut -d " " -f 1 | cut -c1-3)" and resnr "${currentResidueNumber}") or "  >> residuesNotRemodelled.txt
fi
done
fi

# fuse promod3 result and crystal waters and ions, create PDB files to check RMSD
python3 ${PATH_TO_SCRIPTS}/reorganisePDBs.py protein.pdb protein_normalised.pdb missingResidues.txt residuesNotRemodelled.txt

# delete folder if reorganisePDBs.py fails with an index error
if [ -f indexError ]
then
cd ..
rm -rf ${target}
echo "reorganisePDBs.py failed with an index error. The root cause is most likely a subtle error in the gap identification during structure repair. Stopping workflow execution!"
exit 0
fi

# delete folder if modelling is required but not successful
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Stereo-chemical problem in sidechain") && ([ -s missingResidues.txt ] || [ -s missingAtoms.txt ])
then
cd ..
rm -rf ${target}
echo "The crystal structure lacks atoms and could not be modelled successfully. Stopping workflow execution!"
exit 0
fi

# delete folder if a mismatch between the sequence and the structure is suspected
if [ -f atomsInConectStatementsNotFound ]
then
cd ..
rm -rf ${target}
echo "Some atoms in CONECT statements could not be found again after modelling (ScriptError). Stopping workflow execution!"
exit 0
fi

# At this stage, we know whether any residues had to be modelled or atoms in residues were missing
# If both isn't the case, go back to the original structure
# If incomplete residues had to be removed, however, take the modelled structure
if ! [ -s missingResidues.txt ]
then
if ! [ -s missingAtoms.txt ]
then
if ! [ -s residuesNotRemodelled.txt ]
then
echo "CADD: Nothing needed to be modelled. Using the original PDB file!" >> $LOGFILE
echo "Nothing needed to be modelled. Using the original PDB file!"
# final clean up
rm protein.fasta protein.pred protein_normalised.pdb original.pdb remodelled.pdb fused.pdb
mv CONECT_old.txt CONECT.txt
exit 0
fi
fi
fi

# ensure homogeneous numbering in PDB files for RMSD calculation (ignoring gaps)
gmx editconf -f original.pdb -o original2.pdb -resnr 1
mv original2.pdb original.pdb
gmx editconf -f remodelled.pdb -o remodelled2.pdb -resnr 1
mv remodelled2.pdb remodelled.pdb

# check RMSD (in units of Angstrom)
rmsd=$(calculate_rmsd original.pdb remodelled.pdb)
echo "For ${PDBid}, the RMSD between the original crystal structure and the remodelled part (in A) amounts to: ${rmsd}"
if (( $(echo "${rmsd} > 2" | bc -l) ))
then
echo "For ${target}, modelling resulted in a significant deviation from the original crystal structure (more than 2 A). Deleting folder as this protein cannot be handled by the current workflow!"
cd ..
rm -rf ${target}
exit 0
fi
mv fused.pdb protein.pdb
rm original.pdb remodelled.pdb

# final clean up
rm protein.fasta protein.pred protein_normalised.pdb CONECT_old.txt
fi

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload biopython/1.81
module unload tmbed/1.0.0
module unload boost/1.82.0
module unload openmm/7.7.0
module unload ost/2.4
module unload promod/3.3
module unload gromacs/2023.2
module unload rmsd/1.5.1
module unload pdb-tools/2.5.0

export OPENMM_CUDA_COMPILER=""
echo "CADD: Exiting CheckProtein ">> $LOGFILE
