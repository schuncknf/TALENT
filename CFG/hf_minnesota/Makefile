#SHELL = /bin/tcsh
FC = gfortran
OBJECTS= PARS.o gaussquad.o lague.o wavefunction.o m1body.o v2body.o v2test.o meanfields.o iter.o main.o
OPTS = -O3 -ffree-form -ffree-line-length-1024
HFLIBS = -L$../libOffice -llapack -lblas
main : $(OBJECTS)
	$(FC) $(OBJECTS) $(HFLIBS) -o $@
PARS.o : PARS.f90
	$(FC) -c $(OPTS) $<
gaussquad.o : gaussquad.f90
	$(FC) -c $(OPTS) $<
lague.o : lague.f90
	$(FC) -c $(OPTS) $<
wavefunction.o : wavefunction.f90
	$(FC) -c $(OPTS) $<
m1body.o : m1body.f90
	$(FC) -c $(OPTS) $<
v2body.o : v2body.f90
	$(FC) -c $(OPTS) $<
v2test.o : v2test.f90
	$(FC) -c $(OPTS) $<
meanfields.o : meanfields.f90
	$(FC) -c $(OPTS) $<
iter.o : iter.f90
	$(FC) -c $(OPTS) $<
main.o : main.f90
	$(FC) -c $(OPTS) $<
clean:
	rm -f $(OBJECTS)
	rm -f *.mod
