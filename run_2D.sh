#!/bin/bash

L=$1
BETA=$2
M=$3

BETA1000=$(printf "%.0f" $(echo "${BETA} * 1000" | bc))
if [[ ${M} == -* ]]
then
    # negative mass
    M1000=$(printf "m%.0f" $(echo "-(${M}) * 1000" | bc))
else
    # positive mass
    M1000=$(printf "p%.0f" $(echo "${M} * 1000" | bc))
fi

# construct run id
ID=${L}_${BETA1000}_${M1000}

# create the jobs directory
mkdir -p jobs
cd jobs

# create the 2D directory
mkdir -p 2D
cd 2D

# empty the run directory
rm -rf ${ID}
mkdir ${ID}
cd ${ID}

# create directories for data
mkdir -p {ckpoint,topo_charge,vacuum_trace,pion_corr,plots}

# execute the run
../../../schwinger2peD ${L} ${L} 1 ${BETA} 1.0 ${M} 100 20 500
