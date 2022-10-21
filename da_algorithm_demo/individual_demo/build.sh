#!/bin/sh

export CPPFLAGS_CCPL=-DATM_DA
export libname=libda_individual_demo_atm.so
rm *.o *.mod 
make libso

export CPPFLAGS_CCPL=-DOCN_DA
export libname=libda_individual_demo_ocn.so
rm *.o *.mod
make libso

cp libda_individual_demo*.so $WORK_DIR/CCPL_dir/libs/external_procedures/
