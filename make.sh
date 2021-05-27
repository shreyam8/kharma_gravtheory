#!/bin/bash

# Make script for KHARMA
# Used to decide flags and call cmake
# Usage:
# ./make.sh [clean] [cuda] [debug]

# clean: BUILD by re-running cmake, restarting the make process from nothing
#        That is, "./make.sh clean" == "make clean" + "make"
#        Always use 'clean' when switching Release->Debug or OpenMP->CUDA
# cuda:  Build for GPU with CUDA. Must have 'nvcc' in path
# debug: Configure with debug flags: mostly array bounds checks
#        Note most prints/fluid sanity checks are actually *runtime* parameters
# skx:   Compile specifically for Skylake nodes on Stampede2

### Machine-specific configurations ###
if [[ $(hostname) == "toolbox" ]]; then
  HOST=ferrum
else
  HOST=$(hostname -f)
fi

# Kokkos_ARCH options:
# CPUs: WSM, HSW, BDW, SKX, AMDAVX
# ARM: ARMV8, ARMV81, ARMV8_THUNDERX2
# POWER: POWER8, POWER9
# MIC: KNC, KNL
# GPUs: KEPLER35, VOLTA70, TURING75

# TACC resources
# Generally you want latest Intel/IMPI/phdf5 modules,
# On longhorn use gcc7, mvapich2-gdr, and manually-compiled PHDF5
if [[ $HOST == *".frontera.tacc.utexas.edu" ]]; then
  HOST_ARCH="SKX"
fi
if [[ $HOST == *".stampede2.tacc.utexas.edu" ]]; then
  if [[ "$*" == *"skx"* ]]; then
    HOST_ARCH="SKX"
  else
    HOST_ARCH="KNL"
  fi
fi
if [[ $HOST == *".longhorn.tacc.utexas.edu" ]]; then
  HOST_ARCH="POWER9"
  DEVICE_ARCH="VOLTA70"
  PREFIX_PATH="$HOME/libs/hdf5-gcc7-mvapich2"
fi

# Illinois BH cluster
if [[ $HOST == *".astro.illinois.edu" ]]; then
  # When oneAPI works
  #source /opt/intel/oneapi/setvars.sh
  #PREFIX_PATH="$HOME/libs/hdf5-oneapi"

  module load gnu mpich phdf5
  PREFIX_PATH="$MPI_DIR"

  HOST_ARCH="SKX"
fi
# Except BH27/9
if [[ $HOST == "bh29.astro.illinois.edu" ]]; then
  HOST_ARCH="AMDAVX"
fi
if [[ $HOST == "bh27.astro.illinois.edu" ]]; then
  HOST_ARCH="WSM"
fi

# BP's machines
if [[ $HOST == "fermium" ]]; then
  module load mpi
  HOST_ARCH="AMDAVX"
  DEVICE_ARCH="TURING75"
  # My CUDA installs are a bit odd
  EXTRA_FLAGS="-DCUDAToolkit_INCLUDE_DIR=/usr/include/cuda $EXTRA_FLAGS"
fi
if [[ $HOST == "ferrum" ]]; then
  HOST_ARCH="HSW"
  DEVICE_ARCH="INTEL_GEN"
  PREFIX_PATH="$HOME/libs/hdf5-oneapi"
fi
if [[ $HOST == "cinnabar"* ]]; then
  HOST_ARCH="HSW"
  DEVICE_ARCH="KEPLER35"
  PREFIX_PATH="$HOME/libs/hdf5-oneapi"

  if [[ "$*" == *"cuda"* ]]; then
    export NVCC_WRAPPER_DEFAULT_COMPILER=nvc++
    HOST_ARCH="SNB" # Kokkos doesn't detect/set -tp=haswell for nvc++

    #PREFIX_PATH=""
    PREFIX_PATH="$HOME/libs/hdf5-nvhpc"

    # NVHPC CUDA
    #EXTRA_FLAGS="-DCUDAToolkit_INCLUDE_DIR=/opt/nvidia/hpc_sdk/Linux_x86_64/21.3/cuda/include/ $EXTRA_FLAGS"
    # System CUDA
    EXTRA_FLAGS="-DCUDAToolkit_INCLUDE_DIR=/usr/include/cuda $EXTRA_FLAGS"
    # This makes Nvidia chill about old GPUs, but requires a custom nvcc_wrapper
    #export CXXFLAGS="-Wno-deprecated-gpu-targets"
  fi
fi

# If we haven't special-cased already, guess an architecture
# This ends up pretty much optimal on x86 architectures which don't have
# 1. AVX512
# 2. GPUs
if [[ -z "$HOST_ARCH" ]]; then
  if grep GenuineIntel /proc/cpuinfo >/dev/null 2>&1; then
    HOST_ARCH="HSW"
  fi
  if grep AuthenticAMD /proc/cpuinfo >/dev/null 2>&1; then
    HOST_ARCH="AMDAVX"
  fi
fi

# Add some flags only if they're set
if [[ -v HOST_ARCH ]]; then
  EXTRA_FLAGS="-DKokkos_ARCH_${HOST_ARCH}=ON $EXTRA_FLAGS"
fi
if [[ -v DEVICE_ARCH ]]; then
  EXTRA_FLAGS="-DKokkos_ARCH_${DEVICE_ARCH}=ON $EXTRA_FLAGS"
fi
if [[ -v PREFIX_PATH ]]; then
  EXTRA_FLAGS="-DCMAKE_PREFIX_PATH=$PREFIX_PATH $EXTRA_FLAGS"
fi

### Environment ###
if [[ "$(which python)" == *"conda"* ]]; then
  echo "It looks like you have Anaconda loaded."
  echo "Anaconda forces a serial version of HDF5 which makes this compile impossible."
  echo "Deactivate your environment with 'conda deactivate'"
fi
echo "If this is your Anaconda version of Python, deactivate your environment:"
echo "$(which python)"
echo

if [[ "$*" == *"debug"* ]]; then
  TYPE=Debug
else
  TYPE=Release
fi

### Build ###
SCRIPT_DIR=$( dirname "$0" )
cd $SCRIPT_DIR
SCRIPT_DIR=$PWD

# Strongly prefer icc for OpenMP compiles
# I would try clang but it would break all Macs
if which icpc >/dev/null 2>&1; then
  CXX_NATIVE=icpc
  C_NATIVE=icc
  # Avoid warning on nvcc pragmas Intel doesn't like
  export CXXFLAGS="-Wno-unknown-pragmas"
  #export CFLAGS="-qopenmp"
else
  CXX_NATIVE=g++
  C_NATIVE=gcc
fi

# CUDA loop options: MANUAL1D_LOOP > MDRANGE_LOOP, TPTTR_LOOP & TPTTRTVR_LOOP don't compile
# Inner loop must be TVR_INNER_LOOP
# OpenMP loop options for KNL:
# Outer: SIMDFOR_LOOP;MANUAL1D_LOOP;MDRANGE_LOOP;TPTTR_LOOP;TPTVR_LOOP;TPTTRTVR_LOOP
# Inner: SIMDFOR_INNER_LOOP;TVR_INNER_LOOP
if [[ "$*" == *"sycl"* ]]; then
  CXX=icpx
  EXTRA_FLAGS="-DCMAKE_C_COMPILER=icx $EXTRA_FLAGS"
  OUTER_LAYOUT="MANUAL1D_LOOP"
  INNER_LAYOUT="TVR_INNER_LOOP"
  ENABLE_OPENMP="ON"
  ENABLE_CUDA="OFF"
  ENABLE_SYCL="ON"
elif [[ "$*" == *"cuda"* ]]; then
  CXX="$SCRIPT_DIR/external/parthenon/external/Kokkos/bin/nvcc_wrapper"
  OUTER_LAYOUT="MANUAL1D_LOOP"
  INNER_LAYOUT="TVR_INNER_LOOP"
  ENABLE_OPENMP="ON"
  ENABLE_CUDA="ON"
  ENABLE_SYCL="OFF"
else
  CXX="$CXX_NATIVE"
  OUTER_LAYOUT="MDRANGE_LOOP"
  INNER_LAYOUT="SIMDFOR_INNER_LOOP"
  ENABLE_OPENMP="ON"
  ENABLE_CUDA="OFF"
  ENABLE_SYCL="OFF"
fi

# Make build dir. Recall "clean" means "clean and build"
if [[ "$*" == *"clean"* ]]; then
  rm -rf build
fi
mkdir -p build
cd build

if [[ "$*" == *"clean"* ]]; then
#set -x
  cmake ..\
    -DCMAKE_CXX_COMPILER="$CXX" \
    -DCMAKE_INCLUDE_PATH="$PREFIX_PATH" \
    -DCMAKE_BUILD_TYPE=$TYPE \
    -DPAR_LOOP_LAYOUT=$OUTER_LAYOUT \
    -DPAR_LOOP_INNER_LAYOUT=$INNER_LAYOUT \
    -DKokkos_ENABLE_OPENMP=$ENABLE_OPENMP \
    -DKokkos_ENABLE_CUDA=$ENABLE_CUDA \
    -DKokkos_ENABLE_SYCL=$ENABLE_SYCL \
    $EXTRA_FLAGS
#set +x
fi

make -j
cp kharma/kharma.* ..
