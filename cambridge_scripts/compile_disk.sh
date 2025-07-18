#!/bin/bash
#SBATCH --job-name=compile
#SBATCH --ntasks 1
#! Node type
#SBATCH -p genoa
#!How much wallclock time will be required (HH:MM:SS)
#SBATCH --time=00:30:00


cd ../
module purge
module load intel/compiler-intel-llvm/2025.1.0 intel/mpi/2021.15 intel/mkl/2025.1 hdf5/1.14.3/oneapi/openmpi/lmtvkwxw
module load hdf5/1.14.3/oneapi/openmpi/lmtvkwxw

# Actually using g++ for now 
# Need to do further tests for the intel compilers as they typically faster
# python3 configure.py --prob disk --coord cylindrical -omp -hdf5
# python3 configure.py --prob diskplanet_3d_sph --coord spherical_polar --eos adiabatic -omp -hdf5
python3 configure.py --prob diskplanet_2d --coord=cylindrical --eos adiabatic -omp -hdf5
make clean
make

