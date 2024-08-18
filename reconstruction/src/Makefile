# Generated automatically from Makefile.in by configure.
#+
# File: Makefile.in
# Description : Makefile Template for gsim_analysis
#-

# 1. User Specifications                   
PROG = RecoDecays
#PROG = skim_v1
OBJS1 = pi0eta_prob
OBJS2 = ksfwmoments
OBJS3 = BUserInfo
OBJS4 = geninfo
OBJS5 = GammaUserInfo

#KM OBJS = $(PROG).o                                                       
#KM OBJS = $(PROG).o $(OBJS1).o
#KM OBJS = $(PROG).o $(OBJS1).o $(OBJS2).o
OBJS = $(PROG).o $(OBJS1).o $(OBJS2).o $(OBJS3).o $(OBJS4).o $(OBJS5).o
MODULE  = $(PROG).so




#LIBS = -L$(BELLE_RUN_DIR)/lib/so -ltuple

LIBS = -L$(MY_RUN_DIR)/lib/so -L$(BELLE_RUN_DIR)/lib/so -ltuple -lmdst -lkid -lprobutil -lparticle -lhelix -lkfitter -levtvtx -leid -lhamlet -lip -licpv_skim -lbenergy -lnisKsFinder -ltagv -lpntdb -L/bn06/belle/local/pgsql/lib -lcom-svd 


# 2. System Specifications
#    --- Do not change without knowledge

# Compiler Setup with machine dependence

FC = f77
CC = gcc
CXX = c++

DEFS =  -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_POSTGRES=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBTERMCAP=1 -DHAVE_LIBHISTORY=1 -DHAVE_LIBREADLINE=1 -DHAVE_HISTORY=1 -DHAVE_LIBBSD=1 -DHAVE_LIBM=1 -DHAVE_LIBDL=1 -DHAVE_LIBNSL=1 -DHAVE_LIBCRYPT=1 -DHAVE_LIBNSL=1 -DHAVE_LIBDL=1 -DFORTRAN_PPU=1 -DHAVE_LIBCRYPT=1  -DCERNLIB_TYPE

# DEFS =  -DWORDS_BIGENDIAN=1 -DHAVE_LIBCURSES=1 -DHAVE_LIBREADLINE=1 -DHAVE_LIBNSL=1 -DHAVE_LIBSOCKET=1 -DHAVE_LIBDL=1 -DHAVE_LIBPOSIX4=1 -DFORTRAN_PPU=1 -DHAVE_LIBCRYPT=1  -DCERNLIB_TYPE
CPPFLAGS = 
DEPCPPFLAGS = -MM

FFLAGS = -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\"  -fno-second-underscore -fno-automatic -finit-local-zero -fno-emulate-complex
CFLAGS = -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -g -O2
CXXFLAGS = -DHEP_SHORT_NAMES -DBELLE_SHORT_NAMES -DDSTXX_NOINLINE -DBELLE_TARGET_H=\"belle-x86_64-unknown-linux-gnu-g++.h\" -g -fPIC -std=c++0x
SOFLAGS = -shared -Wl,-export-dynamic
LDFLAGS = 

SYSLIB = -lcrypt   -L/belle/local/lib/gcc-lib/i686-pc-linux-gnu/3.0.4 -L/belle/local/lib/gcc-lib/i686-pc-linux-gnu/3.0.4/../../.. -lgfortran -lm -lgcc -lgcc

#SYSLIB = -lcrypt   -L/usr/lib/gcc/x86_64-redhat-linux/3.4.6 -lg2c -lm -lgcc
CLHEPLIB = -lbelleCLHEP

MOTIF_LIBS = -L/usr/X11R6/lib -lXm

POSTGRESINC = /belle/local/include
# Include directories

INCLUDES_C = $(BELLE_TOP_DIR)/include $(NEUROBAYES)/include 
INCLUDES_FORTRAN = $(BELLE_TOP_DIR)/inc

# Dependence description

include $(BELLE_RUN_DIR)/src/config/Makefile.panther

COMPILE_FCPP := $(FC) -c $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(FFLAGS)
COMPILE_FC := $(FC) -c  $(INCLUDES_FORTRAN:%=-I%) $(FFLAGS)
COMPILE_CC := $(CC) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CFLAGS)
COMPILE_CXX := $(CXX) -c  $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(CXXFLAGS)

LINK_FCPP := $(FC)
LINK_FC := $(FC)
LINK_CC := $(CC)
LINK_CXX := $(CXX)

DEPEND_FCPP := $(FC) -M $(DEFS) $(PANTHER_FMACROS) $(INCLUDES_FORTRAN:%=-I%) $(CPPFLAGS) $(filter-out -KPIC -Nq1000,$(CFLAGS))
DEPEND_CC := $(CC) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(filter-out -fPIC -KPIC,$(CFLAGS))
DEPEND_CXX := $(CXX) $(DEPCPPFLAGS) $(DEFS) $(PANTHER_CMACROS) $(INCLUDES_C:%=-I%) $(CPPFLAGS) $(filter-out -fPIC -KPIC,$(CXXFLAGS))

LINUX_G77_BUG = @LINUX_G77_BUG@

%.o:%.c
	$(COMPILE_CC) $<

%.d:%.c
	$(SHELL) -ec '$(DEPEND_CC) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g' > $@'

%.o:%.cc
	$(COMPILE_CXX) $<

%.d:%.cc
	$(SHELL) -ec '$(DEPEND_CXX) $< | sed -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.h/.tdf/g'> $@'


%.o:%.F
	$(COMPILE_FCPP) $<

%.d:%.F
	$(SHELL) -ec '$(DEPEND_FCPP) $< | sed $(LINUX_G77_BUG) -e "s/$*.o[ :]*/$@ &/g" -e 's/\.[12][0-9][0-9][0-9][0-9][0-9][0-9][0-9][a-z]\.inc/.tdf/g'> $@'


# CERNLIB

ifeq "$(CERN)/$(CERN_LEVEL)" "/"
  CERNLIB_LIB_DIR = /belle/cern/2001/lib
else
  CERNLIB_LIB_DIR = $(CERN)/$(CERN_LEVEL)/lib
endif

# Dependencies

all::	$(OBJS)
	$(LINK_CXX) -o $(MODULE) $(SOFLAGS) $(OBJS) $(LIBS) \
	$(CLHEPLIB) $(CERNLIB) $(SYSLIB) 

check:
	./run.csh $(BELLE_LEVEL) > ./out 2>&1

clean::
	rm -f $(OBJS) $(MODULE)  fort.* fpda_pid.* test.hbk panther.dat \
	cdcgeo.bdat* panther?.dat out out? cdchit.paw hit.data hit2.data zhit.data \
	flukaerr.dat last.kumac* track.data

