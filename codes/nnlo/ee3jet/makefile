COMPILER = ifort

ANALYSIS=-thrust-SW

ifeq ($(COMPILER),ifort)
FC = ifort
#FFLAGS = -real-size 128
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

INCLUDE = modules.f90 usermodules.f90 phasespace.f90 kaleu_particles.f90 \
	  statistics.f90 observables.f90 rambo.f90

MODULES = modules.o usermodules.o statistics.o phasespace.o kaleu_particles.o \
	  observables.o rambo.o

OBJF77 = q2g1.o q2g1A50.o q2g1A5i.o \
	 q2g2.o q2g2A60.o q2g2A6i.o q2g2A6st.o q2g2Fax.o q2g2Fcc.o q2g2Fsc.o q2g2Fv.o \
	 q4g0.o \
	 q4g0born.o \
	 q2g2NT.o q2g2A60NT.o \
	 q2g3NT.o q2g3A70NT.o \
	 m_2q2Q1g_tree.o amp_2q2Q1g_tree.o \
	 m_4q1g_tree.o amp_4q1g_tree.o \
	 Lfunctions.o cernlib.o mutil.o \
	 eventshape.o sqmf.o tdhpl.o

USER =  contb.o contr.o contv.o contvv.o contrr.o contrv.o init_processes.o \
	analysis$(ANALYSIS).o \
	born.o \
	real.o \
	virt.o \
	vvirt.o \
	rvirt.o \
	rreal.o \
	cuts.o \
	kaleu-ps.o \
	rambo-ps.o \
	fastjetfortran.o \
	iopee3jet.o \
	iopuser.o \
	$(OBJF77) $(FJCORE)

NNLO = $(MODULES) main.o genprocesses.o boost.o \
       reshuffle.o debug.o random.o regions.o sigma.o alpha.o \
       subtractions.o checklimits.o checkpoles.o ap.o iops.o multi_rambo.o \
       ranlux.o hplog.o regioncheck.o \
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
	cd $(KALEUDIR) ; export FC='$(FC)' ; export USERDIR=$(PWD); export OPT='$(OPT)'; export DEBUG='$(DEBUG)'; make -f makefile_nnlo ; cd -

clean :
	rm nnlo obj/*.o obj/*.a *~ *.mod
