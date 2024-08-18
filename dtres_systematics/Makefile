#
# Makefile for examples
#




MYTARGET_TATAMI = exe/test
MYOBJS_TATAMI = src_test/test.o src_aux/options.o

MYTARGET2_TATAMI = exe/randomize_dtresPars
MYOBJS2_TATAMI = src_systematics/randomize_dtresPars.o src_aux/options.o

MYTARGET3_TATAMI = exe/dtres_syst
MYOBJS3_TATAMI = src_systematics/dtres_syst.o src_aux/options.o

MYTARGET4_TATAMI = exe/signalOnly_LifetimeFit
MYOBJS4_TATAMI = src_test/signalOnly_LifetimeFit.o src_aux/options.o

MYTARGET5_TATAMI = exe/dtres_syst_dalitz
MYOBJS5_TATAMI = src_systematics/dtres_syst_dalitz.o src_aux/options.o

MYTARGET6_TATAMI = exe/belle2_qq_deltaT_trialFit
MYOBJS6_TATAMI = src_test/belle2_qq_deltaT_trialFit.o src_aux/options.o


# Compiler
CXX = g++

# Optimization and Debug flags
OPT = -O2 #-g -v

# Other compiler flags
CXXFLAGS = $(OPT) -Wall -Wextra -std=c++0x

# ROOT flags and libs
ROOT_FLAGS_INC = $(shell root-config --cflags)
ROOT_LIB = $(shell root-config --libs)
ROOFIT_LIB = -lRooFitCore -lRooFit -lHtml -lMinuit -lFoam -lGui
# Link against libMathMore if it exists
ROOFIT_LIB += $(shell if [ -e $(ROOTSYS)/lib/libMathMore.so ]; then echo -lMathMore; fi)

# Boost Library
BOOST_LIB = -lboost_program_options -lboost_serialization

# libRooComplexPDF
LCPDF_INC = -I$(BELLE_FLC)
LCPDF_LIB = -L$(BELLE_FLC)/libComplexPDF -lComplexPDF

# libRooComplexPDF
LRCPDF_INC = -I$(BELLE_FLC) $(LCPDF_INC)
LRCPDF_LIB = -L$(BELLE_FLC)/libRooComplexPDF -lRooComplexPDF $(LCPDF_LIB)

# libRooTatami
LRT_INC = -I$(BELLE_FLC)
LRT_LIB = -L$(BELLE_FLC)/libRooTatami

# All Includes
INCS = -I. $(ROOT_FLAGS_INC) $(LRCPDF_INC) $(LRT_INC)

# All libs 
LIBS = -L. $(ROOT_LIB) $(ROOFIT_LIB) $(BOOST_LIB)

.SUFFIXES: .c .cc

# C++ source files
.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCS) -o $@ $<

# ------------------------------------------------------------------------

all::
ifdef BELLE_LEVEL
	@echo "BELLE_LEVEL "$(BELLE_LEVEL)" detected"
	$(MAKE) $(MYTARGET_TATAMI)
	$(MAKE) $(MYTARGET2_TATAMI)
	$(MAKE) $(MYTARGET3_TATAMI)
	$(MAKE) $(MYTARGET4_TATAMI)
	$(MAKE) $(MYTARGET5_TATAMI)
	$(MAKE) $(MYTARGET6_TATAMI)

else
	@echo "-----------------------------------------------------------------------"
	@echo "Belle software seems not to be installed on this machine. At least the"
	@echo "BELLE_LEVEL environment variable does not exist. Probably you need to"
	@echo "source the bashrc_general or cshrc_general."
	@echo "-----------------------------------------------------------------------"
	@echo -n "Searching for precompiled libRooTatamiPortable.so ... "
ifeq ($(shell if [ -e $(BELLE_FLC)/libRooTatami/libRooTatamiPortable.so ]; then echo -lRooTatamiPortable; fi), -lRooTatamiPortable)
	@echo "OK"
	@echo "Using precompiled libRooTatamiPortable.so to built $(TARGET_TATAMI)"
	@echo "-----------------------------------------------------------------------"
	$(MAKE) $(MYTARGET_TATAMI)
	$(MAKE) $(MYTARGET2_TATAMI)
	$(MAKE) $(MYTARGET3_TATAMI)
	$(MAKE) $(MYTARGET4_TATAMI)
	$(MAKE) $(MYTARGET5_TATAMI)
	$(MAKE) $(MYTARGET6_TATAMI)

else
	@echo "FAILED"
	@echo "Skipping built of $(MYTARGET_TATAMI)"
	@echo "-----------------------------------------------------------------------"
endif
endif



ifndef BELLE_LEVEL
$(MYTARGET_TATAMI): $(MYOBJS_TATAMI)
	$(CXX) -o $@ $(MYOBJS_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable 

$(MYTARGET2_TATAMI): $(MYOBJS2_TATAMI)
	$(CXX) -o $@ $(MYOBJS2_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable

$(MYTARGET3_TATAMI): $(MYOBJS3_TATAMI)
	$(CXX) -o $@ $(MYOBJS3_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable

$(MYTARGET4_TATAMI): $(MYOBJS4_TATAMI)
	$(CXX) -o $@ $(MYOBJS4_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable

$(MYTARGET5_TATAMI): $(MYOBJS5_TATAMI)
	$(CXX) -o $@ $(MYOBJS5_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable

$(MYTARGET6_TATAMI): $(MYOBJS6_TATAMI)
	$(CXX) -o $@ $(MYOBJS6_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatamiPortable

else
$(MYTARGET_TATAMI): $(MYOBJS_TATAMI)
	$(CXX) -o $@ $(MYOBJS_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

$(MYTARGET2_TATAMI): $(MYOBJS2_TATAMI)
	$(CXX) -o $@ $(MYOBJS2_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

$(MYTARGET3_TATAMI): $(MYOBJS3_TATAMI)
	$(CXX) -o $@ $(MYOBJS3_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

$(MYTARGET4_TATAMI): $(MYOBJS4_TATAMI)
	$(CXX) -o $@ $(MYOBJS4_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

$(MYTARGET5_TATAMI): $(MYOBJS5_TATAMI)
	$(CXX) -o $@ $(MYOBJS5_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

$(MYTARGET6_TATAMI): $(MYOBJS6_TATAMI)
	$(CXX) -o $@ $(MYOBJS6_TATAMI) $(LIBS) $(LRT_LIB) -lRooTatami

endif


clean::
	rm -f $(MYOBJS_TATAMI)
	rm -f $(MYOBJS2_TATAMI)
	rm -f $(MYOBJS3_TATAMI)
	rm -f $(MYOBJS4_TATAMI)
	rm -f $(MYOBJS5_TATAMI)
	rm -f $(MYOBJS6_TATAMI)










