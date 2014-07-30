#SHELL = /bin/tcsh
FC = gfortran
OBJECTS= funcmods.o m1body.o v2read.o meanfield.o iter.o init.o test.o 
OPTS = -O3 -ffree-form -ffree-line-length-1024
HFLIBS = -L$../libOffice -llapack -lblas
test : $(OBJECTS)
	$(FC) $(OBJECTS) $(HFLIBS) -o $@
funcmods.o : funcmods.f90
	$(FC) -c $(OPTS) $<
m1body.o : m1body.f90
	$(FC) -c $(OPTS) $<
v2read.o : v2read.f90
	$(FC) -c $(OPTS) $<
meanfield.o : meanfield.f90
	$(FC) -c $(OPTS) $<
init.o : init.f90
	$(FC) -c $(OPTS) $<
iter.o : iter.f90
	$(FC) -c $(OPTS) $<
test.o : test.f90
	$(FC) -c $(OPTS) $<
clean:
	rm -f $(OBJECTS)
	rm -f *.mod
