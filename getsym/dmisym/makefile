# Add file types to the known list of suffixes (file extensions)
.SUFFIXES: .f90 .exe .obj

# The fortran compiler and its options:
PROG = dmisym
FC = ifort
LIBS = /public/home/luohb/source/spglib-1.9.9/spglib-master/lib/libsymspg.a
FFLAGS = -O1 
LFLAGS =

# Define object file names
objs=  spglib_f08.o getsym.o dmisym.o 

.f90.o:
	-$(FC) $(FFLAGS) -c $*.f90

$(PROG).exe: $(objs)
	-$(FC) $(objs) $(LIBS) -o $(PROG).exe
clean:
	rm -f *.o *.mod $(PROG).exe
