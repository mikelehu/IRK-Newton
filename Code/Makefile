#
# make : 20-04-2017 
#
# Description: Makefile for compiling IRK-NEWTON project.
# System: Ubuntu 16.04 lts (64 bits)
# 

# Compiler.
SYS = Linux
CC = /usr/bin/gcc

# Wolfram Mathematica. 
MATHLINKPATH = $(MATHEMATICA_HOME)/SystemFiles/Links/MathLink/DeveloperKit/Linux-x86-64/CompilerAdditions
MPREP = $(MATHLINKPATH)/mprep

# Compiler options.
LIBFLAGS=-L$(MATHLINKPATH)

FLAGS=-I$(MATHLINKPATH) -O2 -Wall -std=c99 -fno-common -mfma
OMP_FLAGS=-fopenmp

# Objects
OBJECTS = Terminal-IRKNEWTON.o  Common-IRKNEWTON.o Problems.o GaussCoefficients.o 
REALOBJECTS = math-IRKNEWTONtm.o math-IRKNEWTON.o Common-IRKNEWTON.o Problems.o GaussCoefficients.o 

.PHONY: clean help 

all : 
	make clean;
	make term-IRKNEWTON;
	make math-IRKNEWTON;

term-IRKNEWTON: $(OBJECTS) 
	${CC} $(OBJECTS) ${OMP_FLAGS} -lm -llapacke -lblas -o Terminal-IRKNEWTON.exe

math-IRKNEWTONtm.c : math-IRKNEWTON.tm
	${MPREP} math-IRKNEWTON.tm -o math-IRKNEWTONtm.c

%.o: %.c def.h 
	${CC}  ${FLAGS} ${OMP_FLAGS} -c -I./  $<

math-IRKNEWTON : $(REALOBJECTS)    
	${CC}  $(REALOBJECTS) $(LIBFLAGS) ${OMP_FLAGS} -lML64i4 -lm -lrt -lstdc++ -luuid -ldl -llapacke -lblas  -o math-IRKNEWTON.exe 

clean :
	rm -rf *.o *~ core *tm.* *.exe 

help:
	@echo "Valid targets;"
	@echo " term-IRKNEWTON"
	@echo " math-IRKNEWTON"




