#!/bin/bash
module purge
module load compiler-intel-llvm/2025.0.4
module load dev-utilities/2025.0.0
module load mpi/2021.14
module load amd/libxc/7.0.0-icx
module load amd/hdf5/2.0.0

export LIB_DIR="lib"
export AMD_OPTS_RELEASE="-O3 -march=native -fma -mprefer-vector-width=512 -qopt-zmm-usage=high"
export AMD_OPTS_DEBUG_C="-g -O0 -traceback"
export AMD_OPTS_DEBUG_F="-g -O0 -traceback -fpe0"
export LIBXC_ROOT="/software/libraries/amd/libxc/7.0.0-icx"
export HDF5_ROOT="/software/libraries/amd/hdf5/2.0.0"

cmake -S . -B build \
    -DCMAKE_INSTALL_PREFIX=$HOME \
    -DCMAKE_C_COMPILER=mpiicx \
    -DCMAKE_Fortran_COMPILER=mpiifx \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_C_FLAGS="$AMD_OPTS_DEBUG_C" \
    -DCMAKE_Fortran_FLAGS="$AMD_OPTS_DEBUG_F -assume byterecl" \
    -DQE_ENABLE_MPI=ON \
    -DQE_ENABLE_OPENMP=ON \
    -DQE_ENABLE_SCALAPACK=ON \
    -DQE_ENABLE_HDF5=ON \
    -DQE_FFTW_VENDOR=Intel_FFTW3 \
    -DQE_BLA_VENDOR=Intel10_64_dyn \
    -DQE_ENABLE_LIBXC=ON \
    -DCMAKE_PREFIX_PATH="${LIBXC_ROOT};${HDF5_ROOT}" \
    -DLIBXC_ROOT=${LIBXC_ROOT} \
    -DLibxc_INCLUDE_DIR=${LIBXC_ROOT}/include \
    -DHDF5_INCLUDE_DIR=${HDF5_ROOT}/include \

# Compile and Staged Install
#make -C build -j 192
#make -C build install DESTDIR=/home/marks/staging_qe_amd
