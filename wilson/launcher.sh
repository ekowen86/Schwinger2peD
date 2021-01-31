#!/bin/bash

export OMP_NUM_THREADS=4

# Simple test script to demonstrate how to use the 2+eD U(1) code

mkdir -p {gauge,data/{data,plaq,creutz,polyakov,rect,top,pion,vacuum,eig}}

# The value of the coupling in the U(1) 2D theory
BETA=3.0

# The value of the coupling in the extra dim
BETAZ=1.0

# The total number of HMC iterations to perform.
HMC_ITER=1000
# The number of HMC iterations for thermalisation.
HMC_THERM=50

# The number of HMC iterations to skip bewteen measurements.
HMC_SKIP=10
# Dump the gauge field every HMC_CHKPT iterations after thermalisation.
HMC_CHKPT=20
# If non-zero, read in the HMC_CHKPT_START gauge field.
HMC_CHKPT_START=0
# HMC time steps in the integration
HMC_NSTEP=30
# HMC trajectory time
HMC_TAU=1.0

# Number of APE smearing hits to perform when measuring topology
APE_ITER=0
# The alpha value in the APE smearing
APE_ALPHA=0.5

# The RNG seed
RNG_SEED=1234

# DYNAMIC (1) or QUENCHED (0)
DYN_QUENCH=1

# LOCKED (1) or DYNAMIC (0)
LOCKED_Z=1

# Dynamic fermion parameters
# Fermion mass
MASS=0.2
# Maximum CG iterations
MAX_CG_ITER=10000
# CG tolerance
CG_EPS=1e-16

# Measuremets: 1 = measure, 0 = no measure
# Polyakov loops
MEAS_PL=1
# Wilson loops and Creutz ratios
MEAS_WL=1
# Pion Correlation function
MEAS_PC=1
# Vacuum trace
MEAS_VT=1

LX=16
LY=16
LZ=1

command="./wilson2peD $BETA $BETAZ $HMC_ITER $HMC_THERM $HMC_SKIP $HMC_CHKPT
    $HMC_CHKPT_START $HMC_NSTEP $HMC_TAU $APE_ITER $APE_ALPHA $RNG_SEED
    $DYN_QUENCH $LOCKED_Z $MASS $MAX_CG_ITER $CG_EPS $MEAS_PC $MEAS_WL $MEAS_VT
    $LX $LY $LZ"

echo $command

$command
