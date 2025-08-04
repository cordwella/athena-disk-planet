#!/bin/bash
#SBATCH --job-name=3d_disk_plan
#SBATCH --ntasks 72
#SBATCH --nodes 1
#! Node type
#SBATCH -p cosmosx
#!How much wallclock time will be required (HH:MM:SS)
#SBATCH --time=12:00:00


module purge
module load intel-oneapi-compilers/2024.1.0/gcc/o3pkwb7f
module load hdf5/1.14.3/oneapi/intel-oneapi-mpi/j3cst467

mpirun ../bin/athena -i athinput.diskplanet_3d