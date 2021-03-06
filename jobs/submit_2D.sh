#!/bin/bash

L=$1
BETA=$2

# mass range
M_START=0.0
M_END=-0.01
M_INC=-0.02

# hmc parameters
N_STEP=30
TAU=1.0

# number of trajectories
N_THERM=200
N_SKIP=20
N_TRAJ=10000

WALL_HOURS=24

for M in $(seq ${M_START} ${M_INC} ${M_END})
do
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

    # enter the 2D directory
    mkdir -p 2D
    cd 2D

    # delete and enter the run directory
    rm -rf ${ID}
    mkdir -p ${ID}
    cd ${ID}

    # create checkpoint directory
    mkdir -p {ckpoint,wf}

    # create a script file with the executable and parameters
    echo "#!/bin/bash" >> run_${ID}.sh
    echo "#$ -l h_rt=${WALL_HOURS}:00:00" >> run_${ID}.sh
    echo "#$ -P qfe" >> run_${ID}.sh
    echo "#$ -N schwinger_${ID}" >> run_${ID}.sh
    echo ../../../schwinger2peD \
        ${L} ${L} 1 ${BETA} 1.0 ${M} \
        ${N_THERM} ${N_SKIP} ${N_TRAJ} \
        ${N_STEP} ${TAU} \
        >> run_${ID}.sh

    # submit to the queue
    echo "submitting job: ${ID}"
	qsub run_${ID}.sh

    # go back to the parent directory
    cd ../..
done
