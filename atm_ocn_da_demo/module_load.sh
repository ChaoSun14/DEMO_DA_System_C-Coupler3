#!/bin/sh


MYPATH=$(readlink -f "${BASH_SOURCE[0]}")
MYPATH=$(dirname "$MYPATH")

module unload mpi/hpcx/2.7.4/gcc-7.3.1
module load compiler/intel/2017.5.239
module load mpi/hpcx/2.7.4/intel-2017.5.239
#module load mpi/intelmpi/2017.4.239
module load mathlib/netcdf/intel/4.4.1
module load mathlib/cdo/intel/1.10.19
module load apps/nco/intel/4.8.1
module load apps/ncl_ncarg/6.3.0
#source /public/software/mpi/intelmpi/2017.4.239/intel64/bin/mpivars.sh intel64
source /opt/hpc/software/compiler/intel/intel-compiler-2017.5.239/bin/compilervars.sh intel64
