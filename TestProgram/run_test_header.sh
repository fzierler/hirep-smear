#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --oversubscribe
#SBATCH -A scw1019
#SBATCH -t 180

module load compiler/intel/2019
module load mpi/intel/2019

