/SHELL = /bin/tcsh

# The C compiler
CC = gcc

# The Fortran compiler
FC = gfortran

# The executable goes here
BINDIR = ./

# Name of executable
EXENAME = test_ho

LIBS = -L$/. -llapack -lblas

OPTS = -O3  -funroll-loops -ffast-math -g -fcheck=all -fbounds-check -ffree-form -ffree-line-length-1024 -fbacktrace -ffpe-trap=invalid -fdefault-double-8 -fdefault-real-8

# Notice that use statements may only invoke modules lying further upp (to the left) 
# in the list
OBJECTS = types.o constants.o sorting.o symmat_mod.o numerical_integration.o Ho_basis.o test_ho.o

#.f90.o :
#	$(FC) -c $(OPTS) $<

$(BINDIR)/$(EXENAME) : $(OBJECTS)
	$(FC) $(OBJECTS) $(LIBS) -o $@

types.o : types.f90
	$(FC) -c $(OPTS) $<

Ho_basis.o : Ho_basis.f90
	$(FC) -c $(OPTS) $<

#Ho_basis.o : Ho_basis.f03
#	$(FC) -c $(OPTS) $<


numerical_integration.o : numerical_integration.f90
	$(FC) -c $(OPTS) $<

test_ho.o : test_ho.f90
	$(FC) -c $(OPTS) $<

symmat_mod.o : symmat_mod.f90
	$(FC) -c $(OPTS) $<

constants.o : constants.f90
	$(FC) -c $(OPTS) $<

sorting.o : sorting.f90
	$(FC) -c $(OPTS) $<

clean:
	rm -f $(OBJECTS)
	rm -f *.mod