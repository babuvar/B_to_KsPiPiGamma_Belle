CERN_ROOT=/sw/belle/cern/2006/
BELLE_LIB_DIR=/sw/belle/belle/b20090127_0910/x86_64-unknown-linux-gnu/opt/lib/

all: calc_observed2true

calc_observed2true: calc_observed2true.o minuit.o
	g++ -g  -o $@ $^ -lm -lgcc -L$(CERN_ROOT)/lib64 -lpacklib -lmathlib -lgfortran -lgcc 

clean::
	rm -f calc_observed2true *.o

%.o: %.c
	gcc -g  -I. -I.. -I../.. -c $< -o $@

%.o: %.cc
	g++ -g  -I. -I.. -I../.. -c $< -o $@

