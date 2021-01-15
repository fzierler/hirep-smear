#!/usr/bin/bash
#SBATCH -n 2
#SBATCH --oversubscribe
#SBATCH -A scw1019
#SBATCH -t 20
#SBATCH -J HRTS # HiRep Test Suite
#SBTACH -o outSED_I.txt
#SBTACH -e errSED_I.txt

module load compiler/intel/2019
module load mpi/intel/2019

