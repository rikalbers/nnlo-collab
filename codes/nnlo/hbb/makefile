COMPILER = ifort

ifeq ($(COMPILER),ifort)
FC = ifort -real-size 64
FFLAGS = -real-size 64
CXX = icpc
CXXFLAGS =
DEBUG = #-g -traceback
OPT = -O3
SAVEFLAG = #-save
STATIC=
LIBS = -cxxlib
endif
ifeq ($(COMPILER),gfortran)
FC = gfortran
FFLAGS = -fdefault-real-8 -fdefault-double-8
CXX = g++
CXXFLAGS = #-stdlib=libc++
DEBUG = -g -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow,underflow,denormal
OPT = -O0
SAVEFLAG = #-fno-automatic
STATIC=
LIBS= -lstdc++

endif

SMEDIR= $(PWD)/smes

PWD = $(shell pwd)
HOME = $(shell dirname $(PWD))

VPATH= ./:../:obj/:$(SMEDIR)/

INCLUDE0 = $(PWD)

FF = $(FC) $(FFLAGS) $(OPT) $(DEBUG) $(STATIC) -I$(INCLUDE0)

CPP = $(CXX) $(OPT) $(CXXFLAGS) -I$(INCLUDE0)

LIBS+= -lkaleu_nnlo

KALEUVER = Kaleu101212
KALEUDIR = $(HOME)/$(KALEUVER)
KALEU = kaleu_modelfile.f90 kaleu_particles.f90 kaleu_userfile.f90 \
	kaleu_random.f #ranlux.f

FJCORE = fjcore.o JadePlugin.o

INCLUDE = modules.f90 phasespace.f90 usermodules.f90 kaleu_particles.f90 \
	  statistics.f90 observables.f90

MODULES = modules.o statistics.o phasespace.o usermodules.o kaleu_particles.o \
	  observables.o 

OBJF77 = cernlib.o

USER =  contb.o contr.o contv.o init_processes.o \
	analysis.o \
	born.o \
	real.o \
	virt.o \
	cuts.o \
	kaleu-ps.o \
	fastjetfortran.o \
	$(OBJF77) $(FJCORE)

NNLO = $(MODULES) main.o genprocesses.o boost.o \
       reshuffle.o debug.o random.o regions.o sigma.o alpha.o \
       checklimits.o subtractions.o ap.o iops.o\
       $(USER)

nnlo : $(NNLO)
	$(FF) $(patsubst %,obj/%,$(NNLO)) -L$(KALEUDIR) $(LIBS) -o $@

%.o : %.f 
	$(FF) -c -o obj/$@ $<

%.o : %.f90 $(INCLUDE)
	$(FF) -c -o obj/$@ $<

%.o : %.cc 
	$(CPP) -c -o obj/$@ $<

kaleu:
	cp $(KALEU) $(PWD)/../$(KALEUVER)/.
	cd $(KALEUDIR) ; export FC=$(FC) ; export USERDIR=$(PWD); export OPT='$(OPT)'; export DEBUG='$(DEBUG)'; make -f makefile_nnlo ; cd -

clean :
	rm nnlo obj/*.o obj/*.a *~ *.mod
