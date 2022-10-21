#!/bin/sh

export CPPFLAGS_CCPL=-DATM_DA
export libname=libda_ensgather_demo_atm.so
rm *.o *.mod 
make libso

export CPPFLAGS_CCPL=-DOCN_DA
export libname=libda_ensgather_demo_ocn.so
rm *.o *.mod
make libso

cp libda_*_demo*.so $WORK_DIR/CCPL_dir/libs/external_procedures/
