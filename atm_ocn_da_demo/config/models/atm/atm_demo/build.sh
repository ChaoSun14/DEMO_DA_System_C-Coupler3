#!/bin/bash

function dump_Macros
{
echo $1
cat > /tmp/makefile.atm << EOF
include $1
all:
	@echo "FC = \${FC}"
	@echo "CC = \${CC}"
	@echo "CXX = \${CXX}"
	@echo "CPP = \${CPP}"
	@echo "FPP = \${FPP}"
	@echo "AR = \${AR}"
	@echo "LD = \${LD}"
	@echo "CFLAGS = \${CFLAGS}"
	@echo "CPPFLAGS = \${CPPFLAGS}"
	@echo "CXXFLAGS = \${CXXFLAGS}"
	@echo "FIXEDFLAGS = \${FIXEDFLAGS}"
	@echo "FREEFLAGS = \${FREEFLAGS}"
	@echo "FFLAGS = \${FFLAGS}"
	@echo "LDFLAGS = \${LDFLAGS}"
	@echo "ULIBS = \${ULIBS}"
	@echo "MOD_SUFFIX = \${MOD_SUFFIX}"
	@echo "NETCDFINC = \${NETCDFINC}"
	@echo "NETCDFLIB = \${NETCDFLIB}"
	@echo "MPIINC = \${MPIINC}"
	@echo "MPILIB = \${MPILIB}"
	@echo "MCPPFLAG = \${MCPPFLAG}"
	@echo "INCLDIR = \${INCLDIR}"
	@echo "SLIBS = \${SLIBS}"
	@echo "CPPDEFS = \${CPPDEFS}"
EOF
    make -f /tmp/makefile.atm >& /tmp/Macros
    atm_file=${RUN_PATH}/src/Makefile.conf
    
    var=$(grep "^NETCDFINC\s*=" /tmp/Macros | sed "s#^NETCDFINC\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^NETCDFINC\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^NETCDFLIB\s*=" /tmp/Macros | sed "s#^NETCDFLIB\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^NETCDFLIB\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^MPIINC\s*=" /tmp/Macros | sed "s#^MPIINC\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^MPIINC\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^MPILIB\s*=" /tmp/Macros | sed "s#^MPILIB\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^MPILIB\s*=\s*\).*#\1 $var#" $atm_file
    
    var=$(grep "^FC\s*=" /tmp/Macros | sed "s#^FC\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^SFC\s*=\s*\).*#\1 $var#" $atm_file
    sed -i "s#\(^SF90\s*=\s*\).*#\1 $var#" $atm_file
    sed -i "s#\(^DM_FC\s*=\s*\).*#\1 $var#" $atm_file
    sed -i "s#\(^DM_F90\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^CC\s*=" /tmp/Macros | sed "s#^CC\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^SCC\s*=\s*\).*#\1 $var#" $atm_file
    sed -i "s#\(^DM_CC\s*=\s*\).*#\1 $var -DMPI2_SUPPORT -DMPI2_THREAD_SUPPORT#" $atm_file

    var=$(grep "^LD\s*=" /tmp/Macros | sed "s#^LD\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^LD\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^CFLAGS\s*=" /tmp/Macros | sed "s#^CFLAGS\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^CFLAGS_LOCAL\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^CPP\s*=" /tmp/Macros | sed "s#^CPP\s*=\(.*\)#\1#")
    sed -i "s#\(^CPP\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^CPPFLAGS\s*=" /tmp/Macros | sed "s#^CPPFLAGS\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^CPPFLAGS_CCPL\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^LDFLAGS\s*=" /tmp/Macros | sed "s#^LDFLAGS\s*=\s*\(.*\)#\1#")
    sed -i "s#\(^LDFLAGS_LOCAL\s*=\s*\).*#\1 $var#" $atm_file

    var=$(grep "^FFLAGS" /tmp/Macros | sed "s# -O[0-3] # #" | sed "s# -i[0-9]* # #" | sed "s# -r[0-9]* # #" | sed "s#^FFLAGS\s*=\(.*\)#\1#")
    sed -i "s#\(^FFLAGS_LOCAL\s*=\s*\).*#\1 $var#" $atm_file

    rm /tmp/makefile.atm
    rm /tmp/Macros
}

# == Get the path of this script ==
MYPATH=$(readlink -f "$0")
MYPATH=$(dirname "$MYPATH")
# =================================
cp $MODEL_DIR/*.f90 ${RUN_PATH}/src/
cp $MODEL_DIR/*.F90 ${RUN_PATH}/src/
cp $MODEL_DIR/Makefile* ${RUN_PATH}/src/

source export_makefile_variables.sh
dump_Macros "$MACFILE"

cd $RUN_PATH/src
make install
exitcode="$?"
if [ "$exitcode" != "0" ]; then
    exit exitcode
fi
cp $RUN_PATH/obj/atm_demo $RUN_PATH/exe/atm_demo
