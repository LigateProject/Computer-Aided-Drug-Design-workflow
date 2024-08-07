#!/bin/bash
set -e

now=$(date "+%d:%m:%Y-%H:%M:%S")
echo "CADD: $now: GROMACS_LiGen_integration_final.sh starting pwd=$PWD" >> $LOGFILE
cp $CADD_INPUT_DIR/$CADD_TARGET/result.csv .
cp $CADD_INPUT_DIR/$CADD_TARGET/docked.mol2 .

# required external software
module use $CADD_SOFTWARE_MODULES
module load python/3.11.1
module load openbabel/3.1.1
module load ambertools/22
module load acpype/2022.7.21
module load gromacs/2023.2

# STaGE itself
module load stage/1.0.0
export GMXLIB=$GMXDATA

target=$(pwd | rev | cut -d "/" -f 1 | rev)

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}
FF=gaff2

PATH_TO_SLURM=${CADD_SUBMISSION_DIR}

echo "Generating topologies for ${target} ..."

echo "Selecting ligands to be used ..."

# should be provided as env variables
resultFile="result.csv"
mol2File="docked.mol2"
numberOfLigands=100 # TODO: turn into user parameter. > actual no. of ligands

# the number of ligands cannot be greater than the number of ligands in the CSV file
if (( ${numberOfLigands} > $(tail -n +2 ${resultFile} | wc -l) ))
then
numberOfLigands=$(tail -n +2 ${resultFile} | wc -l)
echo "Reducing the number of ligands to ${numberOfLigands} because that is the amount of ligands provided after the virtual screen."
fi

# select ligands with highest score
selectedLigands=($(tail -n +2 ${resultFile} | sort -nr -k 3 | cut -d "," -f 1 | head -${numberOfLigands}))

numberOfPoses=1

echo "Generating protein topology ..."
# pdb2gmx changes the atom numbering once again. Let's prepare for that.
if [ -s missingAtoms.txt ]
then
# convert atom numbers to unique descriptions
touch newMissingAtoms.txt
selection=""
numberOfFirstResidue=$(head -1 protein.pdb | cut -c23-27)
searchTerms=$(tail -1 missingAtoms.txt)
for searchTerm in ${searchTerms}; do
findings=$(grep "${searchTerm}" protein.pdb | cut -c$((12-${#searchTerm}))-12)
count=0
for finding in ${findings}; do
count=$((${count} + 1))
if [[ ${finding} == ${searchTerm} ]]
then
break
fi
done

chain=$(grep "${searchTerm}" protein.pdb | head -${count} | tail -1 | cut -c22)
resnr=$(grep "${searchTerm}" protein.pdb | head -${count} | tail -1 | cut -c23-27)
resnr=$(( ${resnr} - ${numberOfFirstResidue} + 1 ))
atomname=$(grep "${searchTerm}" protein.pdb | head -${count} | tail -1 | cut -c14-17)
selection=$(echo ${selection}"(group Protein and chain "${chain}" and resnr "${resnr}" and name "${atomname}") or ")
done

selection=${selection::-4}

echo ${selection} >> newMissingAtoms.txt
mv newMissingAtoms.txt missingAtoms.txt
sed -i missingAtoms.txt -e "s/  )/)/g"
sed -i missingAtoms.txt -e "s/ )/)/g"
fi

file=protein.pdb

# ion residue naming might not meet pdb2gmx's expectations
# TODO: proper fix
sed -i ${file} -e "s/CA   CAS/CA    CA/g"

# check whether the PDB file contains ions or other substances that are not known to the force field used
unknownResidues=("MN" "NI" "CU" "CD" "FE" "CO" "SR")
unknownResidues2=("YCM")
unknownResidueFound=0

for unknownResidue in ${unknownResidues[@]}; do
if grep -q "${unknownResidue}    ${unknownResidue}" ${file} || grep -q "${unknownResidue}   ${unknownResidue}" ${file}
then
cd ..
rm -rf ${target}
echo "An unknown residue (${unknownResidue}) was found in the structure. Stopping workflow execution!"
unknownResidueFound=1
break
fi
done

for unknownResidue in ${unknownResidues2[@]}; do
if grep -q "${unknownResidue}" ${file}
then
cd ..
rm -rf ${target}
echo "An unknown residue (${unknownResidue}) was found in the structure. Stopping workflow execution!"
unknownResidueFound=1
break
fi
done

if (( ${unknownResidueFound} == 1 ))
then
exit 0
fi

# create protein topology
# AMBER99SB-ILDN
# TIP3P
## CONECT records between ions and amino acids are ignored!!!
# identify disulfide bonds
foundDisulfideBond=0
for i in $(seq 1 $(grep "SG" CONECT.txt | wc -l)); do
numberOfSGOccurences=$(grep "SG" CONECT.txt | head -${i} | tail -1 | awk -F"SG" '{print NF-1}')
positionOfFirstSG=$(grep "SG" CONECT.txt | head -${i} | tail -1 | awk -F"SG" '{print$1}')
positionOfFirstSG=${#positionOfFirstSG}
if (( ${numberOfSGOccurences} == 2 )) && (( ${positionOfFirstSG} < 10 ))
then
foundDisulfideBond=1
break
fi
done
if (( ${foundDisulfideBond} == 1 ))
then
gmx pdb2gmx -f $file -renum -ignh -chainsep ter -merge all <<-eof
6
1
eof

# catch long bonds
tailNumber=$(grep -n "Generating topologies for" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Long Bond")
then
echo "The protein has long bonds. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# merge all changes the residue numbering. Update missingResidues.txt and missingAtoms.txt accordingly.
# update residue numbering
firstResidues=$(grep "N3" topol.top | grep "N " | cut -c 18-25 | xargs)

## create lists of chains and the numbers of their initial residues
if $(grep -q "TER" protein.pdb)
then
numberOfChains=$(grep -c "TER" protein.pdb)
chainNames=""
for chain in $(seq 1 ${numberOfChains}); do
lineNumber=$(grep -n "TER" protein.pdb | head -${chain} | tail -1 | cut -d ":" -f 1)
chainNames=$(echo ${chainNames}" "$(head -$((${lineNumber} - 1)) protein.pdb | tail -1 | cut -c22))
done
else
numberOfChains=1
chainNames=$(grep -n "ATOM" protein.pdb | head -1 | cut -d ":" -f 2 | cut -c22)
fi

## loop over files
filenames="missingResidues.txt missingAtoms.txt"
for filename in ${filenames}; do
selection=""
if [ -s ${filename} ]
then
missing=$(tail -1 ${filename})
if $(grep -q "TER" protein.pdb)
then
chainName=""
residueNumber=""
## loop over file content
for word in ${missing}; do
## signal chain name
if [[ $(echo ${selection} | rev | cut -d " " -f 1 | rev) == "chain" ]]
then
chainName=${word}
elif [[ $(echo ${selection} | rev | cut -d " " -f 1 | rev) == "resnr" ]]
then
## formatting
if [[ ${filename} = "missingResidues.txt" ]]
then
word=$(echo ${word} | cut -d ")" -f 1)
fi
## calculate new residue number
residueNumber=$((${word} + $(echo ${firstResidues} | cut -d " " -f $(($(echo ${chainNames} | sed -e "s/ //g" | grep -b -o ${chainName} | cut -d ":" -f 1) + 1))) - 1))
## formatting
if [[ ${filename} = "missingResidues.txt" ]]
then
residueNumber=${residueNumber}")"
fi
word=${residueNumber}
## clear variables after calculating new residue number
chainName=""
residueNumber=""
fi
## append to file content
selection=${selection}" "${word}
done
else
selection=$(tail -1 ${filename})
fi
## remove useless chain identifiers from selection
for chainName in ${chainNames}; do
selection=$(echo ${selection} | sed -e "s/ and chain "${chainName}"//g")
done
## replace files
touch NEW_${filename}
echo ${selection} >> NEW_${filename}
mv NEW_${filename} ${filename}
fi
done

else
gmx pdb2gmx -f $file -renum -ignh -chainsep ter <<-eof
6
1
eof
fi

# catch long bonds
tailNumber=$(grep -n "Generating topologies for" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
tailNumber=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${tailNumber} ))
if $(tail -${tailNumber} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Long Bond")
then
echo "The protein has long bonds. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

ligandCount=0
for lig in ${selectedLigands[@]}; do
echo "Generating ligand topology for ${lig} in lig_${ligandCount} ..."
mkdir lig_${ligandCount}
cd lig_${ligandCount}
# generate topology and .gro files for ligand

## extract poses
mkdir pose_0
python3 ${PATH_TO_SCRIPTS}/extractPose_final.py <<-eof
../${mol2File}
${lig}
pose_0/ligand.mol2
eof
### let obabel adjust the format of the MOL2 file as it needs it
obabel "pose_0/ligand.mol2" -O "pose_0/dummy.mol2"
mv pose_0/dummy.mol2 pose_0/ligand.mol2

mkdir pose_1
python3 ${PATH_TO_SCRIPTS}/extractPose_final.py <<-eof
../${mol2File}
${lig}
pose_1/ligand.mol2
eof
### let obabel adjust the format of the MOL2 file as it needs it
obabel "pose_1/ligand.mol2" -O "pose_1/dummy.mol2"
mv pose_1/dummy.mol2 pose_1/ligand.mol2

## remove LiGen virtual sites
for i in $(seq 0 ${numberOfPoses}); do
cd pose_${i}
if grep -q "EP" ligand.mol2
then
numberOfVirtualSites=$(grep -c "EP" ligand.mol2)
### remove LiGen virtual sites
sed -i '/EP/d' ligand.mol2
### adjust atom number
patternOld="$(head -3 ligand.mol2 | tail -1 | cut -c1-12)"
patternNew="$(( $(head -3 ligand.mol2 | tail -1 | cut -c1-3)-${numberOfVirtualSites} ))""$(head -3 ligand.mol2 | tail -1 | cut -c4-12)"
sed -i ligand.mol2 -e "s/${patternOld}/${patternNew}/g"
### let obabel adjust the format of the MOL2 file as it needs it
obabel "ligand.mol2" -O "dummy.mol2"
mv dummy.mol2 ligand.mol2
fi
cd ..
done

## use MOL2 file of first pose for topology generation
cp pose_0/ligand.mol2 ligand.mol2

# check ligand
## initialise
count=0
countPeptide=0
countDNA=0
countRNA=0
residueCount=0
previousResidueName=""
isPeptideLigand=0
## determine number of atoms
patternNumber=$(grep "TRIPOS" ligand.mol2 | grep -n "@<TRIPOS>ATOM" | cut -d ":" -f 1)
patternNumber=$(( ${patternNumber} + 1 ))
pattern=$(grep -m${patternNumber} "TRIPOS" ligand.mol2 | tail -1)
numberOfAtoms=$(awk '/@<TRIPOS>ATOM$/ {b=NR; next} /'${pattern}'$/ {print NR-b-1; exit}' ligand.mol2)
## determine range for loop over atoms
lineNumber=$(grep -n "@<TRIPOS>ATOM" ligand.mol2 | cut -d ":" -f 1)
lineNumber2=$(( ${lineNumber} + ${numberOfAtoms} ))
## loop over atoms
for i in $(seq -w $(( ${lineNumber} + 1 )) ${lineNumber2}); do
content=$(head -n ${i} ligand.mol2 | tail -n 1 | cut -c59-63)
# remove spaces
content=$(echo ${content})
### check number of residues in ligand
### blind spot: a chain of the same type of residue is considered one residue
### however, pdb2gmx has the same blind spot and will refuse to parameterise a chain containing only one type of residue
if (( ${residueCount} < 2 ))
then
### in the dataset, all hydrogen atoms have a residue name containing ****
### therefore, residue names containing **** must not increase the residue count
if [[ ${content} != ${previousResidueName} ]]
then
residueCount=$(( ${residueCount} + 1 ))
fi
previousResidueName=${content}
fi
### check whether all residues are standard amino acids or nucleobases
### search for the residue name with a space at the end to avoid issues like finding residue "NVA" due to the existence of "NVAL"
### exclude residue names containing numbers
contentCheck=$(echo ${content} | tr -d -c 0-9)
if grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/aminoacids.rtp && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countPeptide=$(( ${countPeptide} + 1 ))
### ugly fix: some small organic compound is abbreviated the same way as an DNA nucleotide
### TODO: find a proper fix
elif grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/dna.rtp && [[ ${content} != "DAN" ]] && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countDNA=$(( ${countDNA} + 1 ))
### ugly fix: some small organic compound is abbreviated the same way as an RNA nucleotide
### TODO: find a proper fix
elif grep -q "${content} " ${GMXDATA}/amber99sb-ildn.ff/rna.rtp && [[ ${content} != "RUN" ]] && (( ${#contentCheck} == 0 ))
then
count=$(( ${count} + 1 ))
countRNA=$(( ${countRNA} + 1 ))
else
break
fi
done

## rename fluorine atoms if unconventional naming is used, e.g. FAA, FAB should be changed to F01, F02
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

if grep -qE -e '^.{8}F' pose_${i}/ligand.mol2
then
countFluorine=1
while IFS= read -r line; do
    if [[ ${line:8:5} =~ F.... ]]; then
        replacement="F$(printf "%-4d" $countFluorine)"
        line="${line:0:8}$replacement${line:13}"
        countFluorine=$((countFluorine+1))
    fi
    echo "$line"
done < pose_${i}/ligand.mol2 > pose_${i}/tmp_mol2
cat pose_${i}/tmp_mol2 > pose_${i}/ligand.mol2
rm pose_${i}/tmp_mol2
fi
done

# (partial) peptide, DNA or RNA ligand
if (( ${countPeptide}>0 )) || (( ${countDNA}>0 )) || (( ${countRNA}>0 ))
then
isPeptideLigand=1
echo "lig_${ligandCount} (${lig}): Found peptide, DNA or RNA residues! The CADD workflow has no special support for these kind of molecules. They are treated like any other small organic molecule (GAFF2 parameterisation)."
fi

# treatment of general ligand
echo "lig_${ligandCount} (${lig}): Assuming a SMALL ORGANIC LIGAND!"
stage=$(which stage.py)
python3 $stage -i "ligand.mol2" -o "ligand_stage" --forcefields $FF

# error catching
## if STaGE fails, try again with different total charge
lineNumberForCheckingStage=$(grep -n "Assuming a SMALL ORGANIC LIGAND!" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
lineNumberForCheckingStage=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${lineNumberForCheckingStage} ))
if $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Error running generator for gaff2") || $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "Error removing OPLS directory")
then
rm -rf *ligand_stage*
python3 $stage -i "ligand.mol2" -o "ligand_stage" --forcefields $FF --alternativeChargeRounding
fi
# if it still fails, give up
lineNumberForCheckingStage=$(grep -n "Assuming a SMALL ORGANIC LIGAND!" ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | tail -1 | cut -d ":" -f 1)
lineNumberForCheckingStage=$(( $(cat ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | wc -l) - ${lineNumberForCheckingStage} ))
if (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error running generator for gaff2") == 2 )) || (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error removing OPLS directory") == 2 ))
then
echo "For lig_${ligandCount} (${lig}), STaGE failed to generate a ligand topology. Deleting folder as this ligand should not be handled by the current workflow."
cd ..
rm -rf lig_${ligandCount}
continue
elif (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error running generator for gaff2") == 1 )) || (( $(tail -${lineNumberForCheckingStage} ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -c "Error removing OPLS directory") == 1 ))
then
echo "For lig_${ligandCount} (${lig}), STaGE worked with the alternative total charge."
else
echo "For lig_${ligandCount} (${lig}), STaGE successfully generated the ligand topology with the total charge from the input MOL2 file."
fi

# sort files
mv ligand_stage_${FF}/ligand_stage.itp ligand_premature.itp
rm posre_ligand_stage.itp

python3 $PATH_TO_SCRIPTS/topologySplitter.py <<-eof
ligand_premature.itp
eof

# write topology summaries to be used for simulations
python3 $PATH_TO_SCRIPTS/writeTopologySummarySingleLigand.py
sed -i topol_amber.top -e "s=topol_=../topol_=g"
sed -i topol_amber.top -e "s=posre=../posre_=g"

# clean-up
rm -rf ligand_stage_${FF} 
rm *.mol2 ligand_stage.gro

# create .gro files for poses and clean up
for i in $(seq 0 ${numberOfPoses}); do

# skip loop iteration if pose was deleted by previous step
if ! [ -d pose_${i} ]
then
continue
fi

obabel "pose_${i}/ligand.mol2" -O "pose_${i}/ligand.gro"
# print .gro file for protein-ligand complex
python3 $PATH_TO_SCRIPTS/printComplexGroFile.py <<-eof
../conf.gro
pose_${i}/ligand.gro
pose_${i}/full.gro
eof
done

# check whether residues very close to the ligand had to be modelled
if [ -s missingResidues.txt ]
then
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
## create .tpr file to produce index file (gmx select requires a .tpr file containing velocities)
gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c pose_${i}/full.gro -r pose_${i}/full.gro -p topol_amber.top -o pose_${i}/check.tpr -po pose_${i}/checkout.mdp -maxwarn 2 || true

### catch potential grompp segfault
if $(tail -1 ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "core dumped")
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), executing gmx grompp resulted in a segmentation fault/core dump. Deleting folder as this complex should not be handled by the current workflow."
rm -rf pose_${i}/
continue
fi

### catch other grompp errors
if [ ! -s pose_${i}/check.tpr ]
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), the TPR file could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -rf pose_${i}/
continue
fi

fi
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) < 2 ));
then
echo "For lig_${ligandCount} (${lig}), the TPR file could be generated for less than two poses. Deleting directory!"
cd ..
rm -rf lig_${ligandCount}
continue
fi

## produce index file (gmx select errors out on composed command => print it to different bash script first)
ligandIdentifier=$(grep "Protein" topol_amber.top | tail -1 | cut -d " " -f 1)
val1=$(grep -n "; Compound" topol_amber.top | cut -d ":" -f 1)
if grep -q "MOL " topol_amber.top
then
val2=$(grep -n "MOL " topol_amber.top | cut -d ":" -f 1)
elif grep -q "${ligandIdentifier} " topol_amber.top
then
val2=$(grep -n "${ligandIdentifier} " topol_amber.top | cut -d ":" -f 1)
fi
diff=$(( ${val2} - ${val1} ))
ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'

neighbours=$(tail -1 missingResidues.txt)
neighbours='"modelledResidues" '${neighbours}'' # promod3 doesn't model H
# keyword chain does not work correctly
if grep -q "Protein_chain" topol_amber.top
then
numberOfChains=$(grep -c "Protein_chain" topol_amber.top)
numberOfChains=$(( ${numberOfChains}/2 ))
for c in $(seq 1 ${numberOfChains}); do
neighbours=$(echo ${neighbours} | sed -e "s/chain $(grep "Protein_chain" topol_amber.top | head -${c} | tail -1 | cut -d "r" -f 2 | cut -d "_" -f 3 | cut -c1)/molecule ${c}/g")
done
else
neighbours=$(echo ${neighbours} | sed -e "s/chain . and //g")
fi
neighbours=\'${neighbours}\'

existingPose="0"

for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i}/check.tpr ]
then
existingPose=${i}
break
fi
done

gmx_command="gmx select -s pose_${existingPose}/check.tpr -on check.ndx -select ${ligand} ${neighbours}"

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh
rm createIndex.sh

## calculate minimum distance between ligand and modelled residues
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
gmx mindist -f pose_${i}/full.gro -s pose_${i}/check.tpr -n check.ndx -od pose_${i}/check.xvg <<-eof
ligand
modelledResidues
eof
distance=$(tail -1 pose_${i}/check.xvg | rev | cut -d " " -f 1 | rev) # distance in nm
distance=$(echo ${distance} 10 | awk '{print $1 * $2}') # distance in A
echo "For pose_${i} of lig_${ligandCount} (${lig}), the minimum distance between the ligand and modelled residues (in A) amounts to: ${distance}"
if (( $(echo "${distance} < 3" | bc -l) ))
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), modelled residues are very close to the bound ligand (less than 3 A). Deleting folder as this structure should not be handled by the current workflow!"
rm -rf pose_${i}
continue
fi

## clean up
rm pose_${i}/checkout.mdp pose_${i}/check.tpr pose_${i}/check.xvg
fi
done
rm check.ndx

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) < 2 ));
then
echo "For lig_${ligandCount} (${lig}), modelled residues are not very close to the bound ligand (less than 3 A) in less than two poses. Deleting directory!"
cd ..
rm -rf lig_${ligandCount}
exit 0
fi

else
echo "No residues needed to be modelled."
fi


# check whether residues very close to the ligand had to be modelled
if [ -s missingAtoms.txt ]
then
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
## create .tpr file to produce index file (gmx select requires a .tpr file containing velocities)
gmx grompp -f ${PATH_TO_SCRIPTS}/mdp/dummy.mdp -c pose_${i}/full.gro -r pose_${i}/full.gro -p topol_amber.top -o pose_${i}/check.tpr -po pose_${i}/checkout.mdp -maxwarn 2 || true

### catch potential grompp segfault
if $(tail -1 ${PATH_TO_SLURM}/slurm-${SLURM_JOB_ID}.out | grep -q "core dumped")
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), executing gmx grompp resulted in a segmentation fault/core dump. Deleting folder as this complex should not be handled by the current workflow."
rm -rf pose_${i}/
continue
fi

### catch other grompp errors
if [ ! -s pose_${i}/check.tpr ]
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), the TPR file could not be generated (grompp error). Deleting folder as this structure should not be handled by the current workflow."
rm -rf pose_${i}/
continue
fi

fi
done

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) < 2 ));
then
echo "For lig_${ligandCount} (${lig}), the TPR file could be generated successfully for less than two poses. Deleting directory!"
cd ..
rm -rf lig_${ligandCount}
exit 0
fi

## produce index file (gmx select errors out on composed command => print it to different bash script first)
ligandIdentifier=$(grep "Protein" topol_amber.top | tail -1 | cut -d " " -f 1)
val1=$(grep -n "; Compound" topol_amber.top | cut -d ":" -f 1)
if grep -q "MOL " topol_amber.top
then
val2=$(grep -n "MOL " topol_amber.top | cut -d ":" -f 1)
elif grep -q "${ligandIdentifier} " topol_amber.top
then
val2=$(grep -n "${ligandIdentifier} " topol_amber.top | cut -d ":" -f 1)
fi
diff=$(( ${val2} - ${val1} ))
ligand='"ligand" molecule '${diff}''
ligand=\'${ligand}\'

neighbours=$(tail -1 missingAtoms.txt)
neighbours='"modelledAtoms" '${neighbours}'' # promod3 doesn't model H
# keyword chain does not work correctly
if grep -q "Protein_chain" topol_amber.top
then
numberOfChains=$(grep -c "Protein_chain" topol_amber.top)
numberOfChains=$(( ${numberOfChains}/2 ))
for c in $(seq 1 ${numberOfChains}); do
neighbours=$(echo ${neighbours} | sed -e "s/chain $(grep "Protein_chain" topol_amber.top | head -${c} | tail -1 | cut -d "P" -f 2 | cut -d "_" -f 3 | cut -c1)/molecule ${c}/g")
done
else
neighbours=$(echo ${neighbours} | sed -e "s/chain . and //g")
fi
neighbours=\'${neighbours}\'

existingPose="0"

for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i}/check.tpr ]
then
existingPose=${i}
break
fi
done

gmx_command="gmx select -s pose_${existingPose}/check.tpr -on check.ndx -select ${ligand} ${neighbours}"

echo "#!/bin/bash" >> createIndex.sh
echo "set -e" >> createIndex.sh
echo "" >> createIndex.sh
echo $gmx_command >> createIndex.sh
bash createIndex.sh
rm createIndex.sh

## calculate minimum distance between ligand and modelled residues
for i in $(seq 0 ${numberOfPoses}); do
if [ -s pose_${i} ]
then
gmx mindist -f pose_${i}/full.gro -s pose_${i}/check.tpr -n check.ndx -od pose_${i}/check.xvg <<-eof
ligand
modelledAtoms
eof
distance=$(tail -1 pose_${i}/check.xvg | rev | cut -d " " -f 1 | rev) # distance in nm
distance=$(echo ${distance} 10 | awk '{print $1 * $2}') # distance in A
echo "For pose_${i} of lig_${ligandCount} (${lig}), the minimum distance between the ligand and modelled atoms (in A) amounts to: ${distance}"
if (( $(echo "${distance} < 3" | bc -l) ))
then
echo "For pose_${i} of lig_${ligandCount} (${lig}), modelled atoms are very close to the bound ligand (less than 3 A). Deleting folder as this protein should not be handled by the current workflow!"
rm -rf pose_${i}
continue
fi

## clean up
rm pose_${i}/checkout.mdp pose_${i}/check.tpr pose_${i}/check.xvg
fi
done
rm check.ndx

# remove directory if none of the poses survives
if (( $(ls -lh | grep "pose_" | wc -l) < 2 ));
then
echo "For lig_${ligandCount} (${lig}), modelled atoms are not very close to the bound ligand (less than 3 A) in less than two poses. Deleting directory!"
cd ..
rm -rf lig_${ligandCount}
exit 0
fi

else
echo "No atoms in residues needed to be modelled."
fi

cd ..
ligandCount=$(( ${ligandCount}+1 ))
done # loop over ligands

# calculate number of ligands for which a topology could be generated
count=$(ls | grep -c "lig")
echo "CADD: ligands created = $count " >> $LOGFILE
if (( ${count} < 2 ))
then
echo "Cannot calculate relative binding free energies with only one ligand. Stopping workflow execution!"
cd ..
rm -rf ${target}
exit 0
fi

# final clean up
rm protein.pdb missingResidues.txt missingAtoms.txt


# check if CONECT records have been ignored (metals are always ignored; only check disulfide bonds)
if [ -s CONECT.txt ]
then

## identify bonds given by CONECT statements
numberOfBonds=$(grep "bond" CONECT.txt | wc -l)
numberOfDisulfideBonds=0
numberOfDisulfideBondsTopology=0
for b in $(seq 1 ${numberOfBonds}); do
bond=$(grep -m${b} "bond" CONECT.txt | tail -1)
atom1=$(echo ${bond} | cut -d " " -f 2)
atom2=$(echo ${bond} | cut -d " " -f 3)

## identify type of atom
### we always have to take the atom type that is listed directly after the first occurrence of the respective atom number in CONECT.txt
element=$(( $(grep -m1 "${atom1}" CONECT.txt | wc -w) - $(grep -m1 "${atom1}" CONECT.txt | awk -F "${atom1}" '{print $2}' | wc -w) + 1 ))
type1=$(grep -m1 "${atom1}" CONECT.txt  | awk '{print $'${element}'}')
element=$(( $(grep -m1 "${atom2}" CONECT.txt | wc -w) - $(grep -m1 "${atom2}" CONECT.txt | awk -F "${atom2}" '{print $2}' | wc -w) + 1 ))
type2=$(grep -m1 "${atom2}" CONECT.txt  | awk '{print $'${element}'}')

## count number of disulfide bonds
if [[ ${type1} == "SG" && ${type2} == "SG" ]]
then
numberOfDisulfideBonds=$((${numberOfDisulfideBonds} + 1))
fi
done # bonds

## identify file containing protein topology
if $(ls | grep -q "_Protein")
then
files=$(ls topol_*.itp)
else
files="$(ls -d */ | head -1)/topol_amber.top"
fi

for file in ${files}; do
## identify all SG in CYS
allSG=""
for s in $(seq 1 $(grep "CYS     SG" ${file} | wc -l)); do
allSG=${allSG}" "$(grep "CYS     SG" ${file} | head -${s} | tail -1 | cut -c1-6 | tr -c -d 0-9)
done

## identify lines where bonds are listed
lineNumber1=$(grep -n "\[ bonds \]" ${file} | cut -d ":" -f 1)
patternNumber=$(grep "\[" ${file} | grep -n "\[ bonds \]" | cut -d ":" -f 1)
patternNumber=$(( ${patternNumber} + 1 ))
pattern=$(grep -m${patternNumber} "\[" ${file} | tail -1 | sed -e "s/\[//g" | sed -e "s/\]//g")
lineNumber2=$(grep -n "\[${pattern}\]" ${file} | cut -d ":" -f 1)

## extract bonds from topology
bondsTopology=$(head -$((${lineNumber2} - 2)) ${file} | tail -$((${lineNumber2} - ${lineNumber1} - 3)))

## search for disulfide bonds
for s in ${allSG}; do
for d in ${allSG}; do
### we have to exclude bonds that have a different first atom whose number just happens to end on the same digits as the atom number of a SG in a disulfide bond
### SG can never be the first atom in a topology such that we can add a leading space to fulfill the above requirement
### similarly, we have to exclude bonds that have a different second atom whose number just happens to begin with the same digits as the atom number of a SG in a disulfide bond
### the bond type is always indicated such that we can add a space after the second atom
if $(echo ${bondsTopology} | grep -q " ${s} ${d} ") && (( ${s} != ${d} ))
then
numberOfDisulfideBondsTopology=$((${numberOfDisulfideBondsTopology} + 1))
fi
done
done

done # files

## PDB CONECT statements list all disulfide bonds twice
numberOfDisulfideBondsTopology=$((${numberOfDisulfideBondsTopology} + ${numberOfDisulfideBondsTopology}))

if (( ${numberOfDisulfideBonds} != ${numberOfDisulfideBondsTopology} ))
then
echo "The number of disulfide bonds in the PDB CONECT statements and in the GROMACS topology is not identical. Check topologies manually before continuing!"
touch incorrect_number_of_disulfide_bond
fi

fi

rm CONECT.txt residuesNotRemodelled.txt

# unload required external software to restore the environment at the beginning of the script
module unload python/3.11.1
module unload openbabel/3.1.1
module unload ambertools/22
module unload acpype/2022.7.21
module unload gromacs/2023.2

export GMXLIB=""
now=$(date "+%d:%m:%Y-%H:%M:%S")
echo "CADD: $now: GROMACS_LiGen_integration_final.sh --exiting  " >> $LOGFILE
