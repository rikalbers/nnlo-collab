FC = gfortran

CUTS = ee
CONSTANTS = ee

RNG = ranlux

KALEUHOME = ..
KALEULIB = kaleu_example

TOYAMPHOME = ./ToyAmp
TOYAMPLIB = toyamp

FFLAGS = 

FILE = montecarlo

all: copy a.out;

%.o : %.f90
	$(FC) $(FFLAGS) -c $*.f90 -o $*.o 

.f.o: 
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

OBJECTS = \
          $(RNG).o \
          ./particles.o \
          ./cuts.o \
          ./rambo.o \
          ./$(FILE).o

a.out: $(OBJECTS) 
	$(FC) $(FFLAGS) -o $(FILE).exe $(OBJECTS) \
	-L$(KALEUHOME) -l$(KALEULIB) \
	-L$(TOYAMPHOME) -l$(TOYAMPLIB)

copy:
	cp cuts_$(CUTS).h90 cuts.h90
	cp constants_$(CONSTANTS).h90 constants.h90

clean:;	rm $(OBJECTS) *.mod *.exe

cuts.o: cuts.h90
montecarlo.o: constants.h90
