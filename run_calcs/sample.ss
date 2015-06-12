#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=gator

#SBATCH --time=0-2:00:00
#SBATCH -e gator.err
#sbatch -o gator.log
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=20

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE


python /directory/to/gator/src/core/master.py -c -f sample.conf -i
