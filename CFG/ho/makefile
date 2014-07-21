#SHELL = /bin/tcsh
FC = gfortran
OBJECTS= pars.o gaussquad.o matrix.o lague.o main.o
OPTS = -O3 -ffree-form -ffree-line-length-1024
HFLIBS = -L$../libOffice -llapack -lblas
hydr : $(OBJECTS)
	$(FC) $(OBJECTS) $(HFLIBS) -o $@
pars.o : pars.f90
	$(FC) -c $(OPTS) $<
gaussquad.o : gaussquad.f90
	$(FC) -c $(OPTS) $<
matrix.o : matrix.f90
	$(FC) -c $(OPTS) $<
lague.o : lague.f90
	$(FC) -c $(OPTS) $<
main.o : main.f90
	$(FC) -c $(OPTS) $<
clean:
	rm -f $(OBJECTS)
	rm -f *.mod
