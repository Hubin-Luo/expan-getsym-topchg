# Add file types to the known list of suffixes (file extensions)
.SUFFIXES: .f90 .exe .obj

# The fortran compiler and its options:
PROG = expan
FC = ifort
LIBS = 
FFLAGS = -O1 
LFLAGS =

# Define object file names
objs=  base.o getpos.o gclst.o  \
       ineqvlt.o correl.o sort.o \
       fit.o sor.o engcal.o main.o 

.f90.o:
	-$(FC) $(FFLAGS) -c $*.f90

$(PROG).exe: $(objs)
	-$(FC) $(objs) $(LIBS) -o $(PROG).exe
clean:
	rm -f *.o *.mod $(PROG).exe
