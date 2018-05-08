#!/bin/bash
#COBALT -n 2 
#COBALT -t 20
#COBALT -A SC_Architectures_2
#COBALT --attrs mcdram=cache:numa=quad

# ^- one note is that lots of Python apps prefer:
#    --attrs mcdram=cache:numa=snc4

# How many ranks to run per node
export PROC_PER_NODE=32

# WRS does this as he loves core files
# when things fail
ulimit -c unlimited

# Activate your environment you're using Intel Python Distribution
#
# WRS recommends cloning your own IPD and if you love performance
# he suggests using the 2018 Beta or asking about ALCF Python
source activate idp

# Setting your number of OpenMP threads depends on your application
# one size does not fit all.
export OMP_NUM_THREADS=1

# This should pick up your VirtualEnv or Conda python
#BINARY="$(which python)"
BINARY="/home/fcurtis/anaconda2/envs/idp/bin/python"

GATOR_MASTER=../src/GAtor_master.py
CONF_FILE="ui.conf"
echo "##########################################"

echo "Running: ${BINARY} ${PYTHON_SCRIPT}"

echo "##########################################"

# Since Intel Python uses Intel MPI, we need to use the ABI compatibility
# See http://docs.cray.com/books/S-2544-704/S-2544-704.pdf for details
module swap cray-mpich cray-mpich-abi
export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:$LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

# What are all these flags?
#   -b :: don't copy the binary
# --cc :: cpu binding mode
#   -d :: cores to bind a rank to
#   -j :: hyperthreading
#   -n :: how many ranks for the job
#   -N :: how many ranks per node
aprun -b -cc depth -d 1 -j 1 -n $(($COBALT_PARTSIZE * $PROC_PER_NODE)) -N $PROC_PER_NODE ${BINARY} ${GATOR_MASTER} ${CONF_FILE}
