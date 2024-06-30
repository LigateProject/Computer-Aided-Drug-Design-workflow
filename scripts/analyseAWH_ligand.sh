#!/bin/bash
set -e

PATH_TO_SCRIPTS=${CADD_SCRIPTS_DIR}

module use $CADD_SOFTWARE_MODULES
module load gromacs/2023.2

read_in_awh () {
    local g=$(grep $1 awh_pmf_ligand_t$2.xvg)
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

awh_nstout=100

# analyse
AVERAGE=0.0
STDDEV=0.0
for pose in pose*; do
cd ${pose}
if (( $(pwd | rev | cut -d "_" -f 1 | rev) == 0 ))
then
LAMBDA0=" 0.0000"
LAMBDA1="47.0000"
else
LAMBDA0="47.0000"
LAMBDA1=" 0.0000"
fi
lineTimeLigand=$(($(grep -n "Step           Time" production_ligand.log | tail -1 | cut -d ":" -f 1)+1 ))
timeLigand=$(head -${lineTimeLigand} production_ligand.log | tail -1 | awk '{print $2}' | cut -d "." -f 1)
timeLigand=$(( ${timeLigand} - $(( ${timeLigand}%${awh_nstout} )) ))
gmx awh -f production_ligand.edr -s production_ligand.tpr -o awh_pmf_ligand.xvg -more -b ${timeLigand} -e ${timeLigand}
DELTAG0=$(read_in_awh ${LAMBDA0} ${timeLigand})
DELTAG1=$(read_in_awh ${LAMBDA1} ${timeLigand})
DELTAG=$(echo "${DELTAG1} ${DELTAG0}" | awk '{print $1-$2}')
AVERAGE=$(echo "${AVERAGE} ${DELTAG}" | awk '{print $1+$2}')
cd ..
done
AVERAGE=$(print_average ${AVERAGE} $(ls | grep -c "pose"))
for pose in pose*; do
cd ${pose}
if (( $(pwd | rev | cut -d "_" -f 1 | rev) == 0 ))
then
LAMBDA0=" 0.0000"
LAMBDA1="47.0000"
else
LAMBDA0="47.0000"
LAMBDA1=" 0.0000"
fi
lineTimeLigand=$(($(grep -n "Step           Time" production_ligand.log | tail -1 | cut -d ":" -f 1)+1 ))
timeLigand=$(head -${lineTimeLigand} production_ligand.log | tail -1 | awk '{print $2}' | cut -d "." -f 1)
timeLigand=$(( ${timeLigand} - $(( ${timeLigand}%${awh_nstout} )) ))
DELTAG0=$(read_in_awh ${LAMBDA0} ${timeLigand})
DELTAG1=$(read_in_awh ${LAMBDA1} ${timeLigand})
DELTAG=$(echo "${DELTAG1} ${DELTAG0}" | awk '{print $1-$2}')
STDDEV=$(add_to_stddev ${STDDEV} ${AVERAGE} ${DELTAG})
STDDEV=$(print_stddev ${STDDEV} $(echo "$(ls .. | grep -c 'pose_*') $(ls .. | grep -c 'pose_*')" | awk '{print $1*$2}'))
rm awh_pmf_ligand*xvg
cd ..
done
echo "${timeLigand} ${AVERAGE} ${STDDEV}" >> freeEnergySummaryLigand.txt

# Check for convergence if we have three or more data points
if (( $(cat freeEnergySummaryLigand.txt | wc -l)>2 ))
then
convergence=$(python3 ${PATH_TO_SCRIPTS}/is_converged.py freeEnergySummaryLigand.txt --max_error 4) # TODO: turn into user parameter
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
fi
