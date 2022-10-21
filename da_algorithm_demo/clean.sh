#!/bin/sh

export ROOT_DIR=/public/home/zhangxz/dafcc_test/sample_da_dafcc1/
export LIB_DIR=$ROOT_DIR/da_algorithm_demo

export CASE_LOCAL_INCL="-check all"
echo " "
echo " "
echo "++++++++++++++++"
echo "Individual da demo clean"
echo "++++++++++++++++"
echo " "
echo " "
cd $LIB_DIR/individual_demo
cp ../Makefile.conf .
make clean
ls
echo " "
echo " "
echo "++++++++++++++++"
echo "Ensmean da demo clean"
echo "++++++++++++++++"
echo " "
echo " "
cd $LIB_DIR/ensmean_demo
cp ../Makefile.conf .
make clean
ls
echo " "
echo " "
echo "++++++++++++++++"
echo "Ensgather da demo clean"
echo "++++++++++++++++"
cd $LIB_DIR/ensgather_demo
cp ../Makefile.conf .
make clean
ls
echo " "
echo " "




