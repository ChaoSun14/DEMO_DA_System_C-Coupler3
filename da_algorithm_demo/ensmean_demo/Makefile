SHELL=/bin/sh

#==============================================================================
#
# da_algorithm_demo Makefile
#
# <Usage>
#
#   0) Modify the include link to either use Makefile.conf.linux
#      or makefile.conf.sgi for compilation on the linux or sgi
#
#   1) To make a executable file, type
#         > make  or  > make all
#
#   2) To make a executable file with debug options, type
#         > make debug
#
#   3) To copy the load module to installing directory, type
#         > make install
#      . Specify the directory to a variable 'INSTALL_DIR' below.
#
#   4) To make dynamic link library, type
#         > make libso
#
#   4) To clean up files created by make, type
#         > make clean
#
# Created by Sun Chao in Sep,2021
#==============================================================================

#-----------------------------------------------------------------------------
#                          -- Parent make (calls child make) --
#-----------------------------------------------------------------------------

# -----------------------------------------------------------
# Default configuration, possibily redefined in Makefile.conf
# -----------------------------------------------------------

ARCH	 = `uname -s`
SED      = sed
DASPERL  = /usr/bin/perl

# -------------
# General Rules
# -------------

CP              = /bin/cp -p
RM              = /bin/rm -f
MKDIR           = /bin/mkdir -p
AR              = ar cq
PROTEX          = protex -f # -l
ProTexMake      = protex -S # -l
LATEX           = pdflatex
DVIPS           = dvips

# Preprocessing
# -------------
_DDEBUG =
_D      = $(_DDEBUG)

# ------------------------
# Default System Libraries
# ------------------------
LIBmpi          = -lmpi
LIBsys          = 


#------------
# Include machine dependent compile & load options
#------------
  MAKE_CONF = Makefile.conf
include $(MAKE_CONF)


# -------------
# This makefile
# -------------

  MAKE_FILE = Makefile


# -----------
# Load module
# -----------

  EXE_FILE = da_algorithm_demo
  #LIB_SO = libda_ensmean_demo_ocn.so
  LIB_SO = $(libname)

# --------------------
# Installing directory
# --------------------

  INSTALL_DIR = ../obj


# --------
# Log file
# --------

  LOG_FILE = log.make.$(EXE_FILE)


# ---------------
# Call child make
# ---------------

"" :
	@$(MAKE) -f $(MAKE_FILE) all


# ------------
# Make install
# ------------

install:
	@echo
	@echo '==== INSTALL ================================================='
	@if [ -e $(INSTALL_DIR) ]; then \
	  if [ ! -d $(INSTALL_DIR) ]; then \
	    echo '### Fail to create installing directory ###' ;\
	    echo '### Stop the installation               ###' ;\
	    exit ;\
	  fi ;\
	else \
	  echo "	mkdir -p $(INSTALL_DIR)" ;\
	  mkdir -p $(INSTALL_DIR) ;\
	fi
	@$(MAKE) -f $(MAKE_FILE) all
	mv *.o *.mod log.* $(INSTALL_DIR)
	cp $(EXE_FILE) $(INSTALL_DIR)
	@cd $(INSTALL_DIR) ; ls -l `pwd`/$(EXE_FILE)


# ----------
# Make clean
# ----------

clean:
	@echo
	@echo '==== CLEAN ==================================================='
	- $(RM) $(EXE_FILE) *.o *.mod *.MOD *.lst *.a *.x *__genmod* *.so
	- $(RM) log.make.$(EXE_FILE)


#-----------------------------------------------------------------------------
#                          -- Child make --
#-----------------------------------------------------------------------------

# ------------
# Source files
# ------------

  SRCSF90C = \
	da_algorithm_demo.F90 \
	ran_mod.F90 \
	decomp_init_mod.F90 \
	grid_init_mod.F90 \
	model_setting_mod.F90 \
	parse_namelist_mod.F90 \
	spmd_init_mod.F90 \
	variable_mod.F90

  SRCSF90C_SO = \
	da_ensmean_demo_mod.F90 \
	da_ensmean_demo_ccpl_coupling_mod.F90 \
	ran_mod.F90 \
	decomp_init_mod.F90 \
	grid_init_mod.F90 \
	model_setting_mod.F90 \
	parse_namelist_mod.F90 \
	spmd_init_mod.F90 \
	variable_mod.F90



  SRCSF77 =

  SRCSC =

  SRCS = $(SRCSF90C) $(SRCSF77) $(SRCSC) 

# ------------
# Object files
# ------------

  SRCSF90	    = ${SRCSF90C:.F90=.f90}
  SRCSF90_SO    = ${SRCSF90C_SO:.F90=.f90}

  OBJS 		    = ${SRCSF90:.f90=.o} ${SRCSF77:.f=.o} ${SRCSC:.c=.o}
  OBJS_SO       = ${SRCSF90_SO:.f90=.o} ${SRCSF77:.f=.o} ${SRCSC:.c=.o}


# ------------
# Dependencies
# ------------
  MAKE_DEPEND = Makefile.dependency
include $(MAKE_DEPEND)

# ----

$(EXE_FILE) :  $(OBJS) 
	$(LD) $(FFLAGS) $(CFLAGS) $(OBJS) $(LIBS) $(LDFLAGS) -fPIC -Wl,--export-dynamic -o $@ 

$(LIB_SO) : $(OBJS_SO) 
	$(LD) $(OBJS_SO) $(LIBS) $(LDFLAGS) -shared -fPIC -o $@
# ------------------------
# Call compiler and linker
# ------------------------

all :
	@$(MAKE) -f $(MAKE_FILE) "COMP_MODE=$@" check_mode
	@echo '==== COMPILE ================================================='
	@$(MAKE) -f $(MAKE_FILE) \
        "FFLAGS=$(FFLAGS)" \
        "CFLAGS=$(CFLAGS)" \
        $(OBJS)
	@echo
	@echo '==== LINK ===================================================='
	@$(MAKE) -f $(MAKE_FILE) \
		"LIBS=$(LIBS)" "LDFLAGS=$(LDFLAGS)" \
		$(EXE_FILE)

debug :
	@$(MAKE) -f $(MAKE_FILE) "COMP_MODE=$@" check_mode
	@echo
	@echo '========================================================'
	@$(MAKE) -f $(MAKE_FILE) \
		"LIBS=$(LIBS_D)" "LDFLAGS=$(LDFLAGS_D)" \
		$(LIB_SO) $(EXE_FILE)

check_mode :
	@if [ -e $(LOG_FILE) ]; then \
	  if [ '$(COMP_MODE)' != `head -n 1 $(LOG_FILE)` ]; then \
	    echo ;\
	    echo "### COMPILE MODE WAS CHANGED ###" ;\
	    $(MAKE) -f $(MAKE_FILE) clean ;\
	  fi ;\
	else \
	  echo ;\
	  echo "### NO LOG FILE ###" ;\
	  $(MAKE) -f $(MAKE_FILE) clean ;\
	fi
	@echo $(COMP_MODE) > $(LOG_FILE)

libso :
	@$(MAKE) -f $(MAKE_FILE) "COMP_MODE=$@" check_mode
	@echo '==== COMPILE ================================================='
	@$(MAKE) -f $(MAKE_FILE) \
        "FFLAGS=$(FFLAGS)" \
        "CFLAGS=$(CFLAGS)" \
        $(OBJS_SO)
	@echo
	@echo '==== LINK ===================================================='
	@$(MAKE) -f $(MAKE_FILE) \
		"LIBS=$(LIBS)" "LDFLAGS=$(LDFLAGS)" \
		$(LIB_SO)

help:
	@ echo "Available targets:"
	@ echo "  make             creates executable        "
	@ echo "  make debug       created exec for debugging purposes"
	@ echo "  make install     creates exec & places it in bin"
	@ echo "  make libso       creates dynamic link library"

