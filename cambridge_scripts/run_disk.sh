#!/bin/bash
#SBATCH --job-name=disk
#SBATCH --ntasks 1
#! One node only, using just OpenMP rather than MPI
#SBATCH --nodes 1
#! Max is 186 on genoa
#SBATCH --cpus-per-task 24 
#! Node type
#SBATCH -p genoa
#!How much wallclock time will be required (HH:MM:SS)
#SBATCH --time=12:00:00

# For simplicity use the 
module purge
module load intel/compiler-intel-llvm/2025.1.0 intel/mpi/2021.15 intel/mkl/2025.1 hdf5/1.14.3/oneapi/openmpi/lmtvkwxw
module load hdf5/1.14.3/oneapi/openmpi/lmtvkwxw # hdf5/1.14.3/oneapi/intel-oneapi-mpi/lxttf5co
export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi2.so
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASKS

../bin/athena -i athinput.disk_cyl -t 11:50:00

