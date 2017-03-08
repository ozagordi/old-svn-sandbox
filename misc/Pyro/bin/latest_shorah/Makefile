SHELL = /bin/bash

# common variables
#HOSTNAME = echo hostname
hostname := $(shell hostname)
whoami := $(shell whoami)
host-type := $(shell uname)
CCFLAGS = -pedantic -Wall -Wextra -O2 -ffast-math -funroll-loops# -m32 -DDEBUG
WFLAGS =
VPATH = ./
OLIBS = 
CPP = g++


# EDIT THESE LINES TO INCLUDE GSL
CXLIBS = -lgsl -lgslcblas
CFLAGS = $(CCFLAGS) -I/opt/local/include
XLIBS  = -L/opt/local/lib $(CXLIBS)
CPPFLAGS = $(CFLAGS)

#SRC_1 = dpm_src/dmm_sampler.cpp
OBJS_1 = dpm_src/dpm_sampler.o
EXE_1 = diri_sampler


%.o : %.cpp %.h data_structures.h Makefile # only used for C program diri_sampler
	@echo ''
	@echo '*********************************'
	@echo '   making object: $@ '
	@echo '*********************************'
	$(CPP) $(CFLAGS) $(WFLAGS) $(XLIBS) -c $< -o $@

all: $(EXE_1)

$(EXE_1): $(OBJS_1) Makefile #diri_sampler
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(XLIBS) $(OBJS_1) -O2 -o $(EXE_1)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

.PHONY : clean

clean:
	rm -rf $(OBJS_1) $(EXE_1) $(EXE_2) $(EXE_3) $(EXE_2).dSYM $(EXE_3).dSYM *pyc ./*pyc
