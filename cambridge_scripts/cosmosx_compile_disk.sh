#!/bin/bash
#SBATCH --job-name=compile
#SBATCH --ntasks 1
#! Node type
#SBATCH -p cosmosx
#!How much wallclock time will be required (HH:MM:SS)
#SBATCH --time=00:30:00


cd ../
module purge
module load intel-oneapi-compilers/2024.1.0/gcc/o3pkwb7f
module load hdf5/1.14.3/oneapi/intel-oneapi-mpi/j3cst467

# python3 configure.py --prob=diskplanet_3d_sph --coord=spherical_polar -hdf5 -mpi --cxx icc --eos adiabatic --flux roe
python3 configure.py --prob=diskplanet_2d --coord=cylindrical -hdf5  -mpi --cxx icc --eos adiabatic --flux roe
make clean
make

