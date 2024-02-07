#!/bin/bash
set -e

PATH_TO_SCRIPTS=CADD/scripts

module load gromacs/2023.2 gromacs=gmx

LAMBDA0=" 0.0000"
LAMBDA1="11.0000"

read_in_awh () {
    local g=$(grep $1 awh_pmf_t$2.xvg)
    local g=$(echo ${g} | awk '{print $2}')
    echo ${g}
}

add_to_stddev () {
    local s=$1
    DIFF=$(echo "$2 $3" | awk '{print $1-$2}')
    DIFF=$(echo "${DIFF} ${DIFF}" | awk '{print $1*$2}')
    local s=$(echo "${s} ${DIFF}" | awk '{print $1+$2}')
    echo ${s}
}

print_average () {
    local a=$1
    local a=$(echo "${a} $2" | awk '{print $1/$2}')
    echo ${a}
}

print_stddev () {
    local s=$1
    local s=$(echo "${s} $2" | awk '{print $1/$2}')
    local s=$(echo "sqrt(${s})" | bc)
    echo ${s}
}

# analyse
AVERAGE=0.0
STDDEV=0.0
for pose in pose*; do
cd ${pose}
rm readyForAnalysisLigand
lineTimeLigand=$(($(grep -n "Step           Time" production_ligand.log | tail -1 | cut -d ":" -f 1)+1 ))
timeLigand=$(head -${lineTimeLigand} production_ligand.log | tail -1 | awk '{print $2}')
gmx awh -f production_ligand.edr -s production_ligand.tpr -o awh_pmf_ligand.xvg -more -b ${timeLigand} -e ${timeLigand}
DELTAG0=$(read_in_awh ${LAMBDA0} ${timeLigand})
DELTAG1=$(read_in_awh ${LAMBDA1} ${timeLigand})
DELTAG=$(echo "${DELTAG1} ${DELTAG0}" | awk '{print $1-$2}')
AVERAGE=$(echo "${AVERAGE} ${DELTAG}" | awk '{print $1+$2}')
cd ..
done
AVERAGE=$(print_average ${AVERAGE} $(ls | grep -c "pose_*"))
for pose in pose*; do
cd ${pose}
DELTAG0=$(read_in_awh ${LAMBDA0} ${timeLigand})
DELTAG1=$(read_in_awh ${LAMBDA1} ${timeLigand})
DELTAG=$(echo "${DELTAG1} ${DELTAG0}" | awk '{print $1-$2}')
STDDEV=$(add_to_stddev ${STDDEV} ${AVERAGE} ${DELTAG})
STDDEV=$(print_stddev ${STDDEV} $(echo "${REPLICATES} ${REPLICATES}" | awk '{print $1*$2}'))
rm awh_pmf_ligand*xvg
cd ..
done
echo "${timeLigand} ${AVERAGE} ${STDDEV}" >> freeEnergySummaryLigand.txt

# Check for convergence if we have three or more data points
if (( $(cat freeEnergySummaryLigand.txt | wc -l)>2 ))
then
convergence=$(python3 is_converged.py freeEnergySummaryLigand.txt --max_error 4) # TODO: user-defined parameter
if [[ ${convergence} == "YES" ]]
then
touch isConvergedLigand
for pose in pose*; do
cd ${pose}
rm pid_ligand.txt
cd ..
done
exit 0
fi

# re-start simulation if not converged
for pose in pose*; do
cd ${pose}
rm pid_ligand.txt
bash ${PATH_TO_SCRIPTS}/monitorAWH_ligand.sh
echo $! >> pid_ligand.txt
cd ..
done
