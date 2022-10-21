#!/bin/sh

export ROOT_DIR=/public/home/zhangxz/dafcc_test/sample_da_dafcc1/
export LIB_DIR=$ROOT_DIR/da_algorithm_demo
export WORK_DIR=$ROOT_DIR/atm_ocn_da_demo
export CCPL_LIB=$WORK_DIR/run/libs/c-coupler/cd3c1c1d44f9340d6b8e9f6b6e1fd6d3c2dce75a
export NETCDFPATH=/public/software/mathlib/libs-intel/netcdf/4.4.1
export MPIPATH=/opt/hpc/software/mpi/hpcx/v2.7.4/intel-2017.5.239/
#export MPIPATH=/public/software/mpi/intelmpi/2017.4.239


export ATM_DA=True #False
export OCN_DA=True #False

#export CASE_LOCAL_INCL="-check all"
echo " "
echo " "
echo "++++++++++++++++"
echo "Individual da demo build"
echo "++++++++++++++++"
echo " "
echo " "
cd $LIB_DIR/individual_demo
cp -rf ../Makefile.conf .
./build.sh
echo " "
echo " "
echo "++++++++++++++++"
echo "Ensmean da demo build"
echo "++++++++++++++++"
echo " "
echo " "
cd $LIB_DIR/ensmean_demo
cp -rf ../Makefile.conf .
./build.sh
echo " "
echo " "
echo "++++++++++++++++"
echo "Ensgather da demo build"
echo "++++++++++++++++"
cd $LIB_DIR/ensgather_demo
cp -rf ../Makefile.conf .
./build.sh
echo " "
echo " "




