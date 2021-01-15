#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --oversubscribe
#SBATCH -A scw1019
#SBATCH -t 20
#SBATCH -J HRTS # HiRep Test Suite
#SBATCH -o outSED_I.txt
#SBATCH -e errSED_I.txt

module load compiler/intel/2019
module load mpi/intel/2019

