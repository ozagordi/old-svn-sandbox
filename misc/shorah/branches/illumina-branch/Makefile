SHELL = /bin/bash

# common variables
#HOSTNAME = echo hostname
hostname := $(shell hostname)
whoami := $(shell whoami)
host-type := $(shell uname)
CC = gcc
CPP = g++
CCFLAGS = -Wall -W -Wextra -Winline -ggdb  -O2 -ffast-math  # -funroll-loops# -m32
WFLAGS =
VPATH = ./
OLIBS = 

# EDIT THESE LINES TO INCLUDE GSL
CXLIBS = -lgsl -lgslcblas
CFLAGS = $(CCFLAGS) -I/opt/local/include
XLIBS  =  $(CXLIBS) -L/opt/local/lib	



OBJS_1 = dpm_src/dpm_sampler.o
EXE_1 = diri_sampler

SRC_2 = contain_src/contain.cc
EXE_2 = contain
FLAGS_2 = -g -O3# -Wall

#OBJS_3 = freqEst_src/freqEst.o
SRC_3 = freqEst_src/freqEst.cc
EXE_3 = freqEst
FLAGS_3 = -g -O3 -Wall -ffast-math


%.o : %.c %.h data_structures.h Makefile # only used for C program diri_sampler
	@echo ''
	@echo '*********************************'
	@echo '   making object: $@ '
	@echo '*********************************'
	$(CC) $(CFLAGS) $(WFLAGS) $(OLIBS) -c $< -o $@

all: $(EXE_1) $(EXE_2) $(EXE_3)

$(EXE_1): $(OBJS_1) Makefile #diri_sampler
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CC)  $(OBJS_1) $(XLIBS) -o $(EXE_1)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

$(EXE_2): Makefile #contain
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(FLAGS_2) $(SRC_2) -o $(EXE_2)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

$(EXE_3): Makefile #freqEst
	@echo ''
	@echo '*********************************'
	@echo ' making executable: $@ '
	@echo '*********************************'
	$(CPP) $(FLAGS_3) $(SRC_3) -o $(EXE_3)

	@echo '*******************'
	@echo 'compiled for $(host-type)'

.PHONY : clean

clean:
	rm -f $(OBJS_1) $(EXE_1) $(EXE_2) $(EXE_3)

