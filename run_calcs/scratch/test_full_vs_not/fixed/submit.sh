#!/bin/bash
#SBATCH --qos=normal
#SBATCH --job-name=fixed
#SBATCH --time=0-06:00:00
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=20
#SBATCH --export =ALL


echo Start Job


export OMP_NUM_THREADS=1

mpirun /home/fcurtis/fhi-aims.150205/bin/aims.150205.scalapack.mpi.x >aims.out

echo End Job

