FC = gfortran $(FFLAGS) 

p1 = run_HF

FFLAGS =  -O3 -fopenmp 
TFLAGS =  -g  -fbounds-check

LIBS =  -L/user/local/lib/ -llapack -lblas

FERN: ME_minn.o Hartree_Fock.o main.o
	${FC} ME_minn.o Hartree_Fock.o main.o  -o ${p1}  ${LIBS} 

FAB: fabrizio.f90
	${FC} fabrizio.o -o ${p3} ${LIBS}

ME_minn.o: ME_minn.f90 
	${FC} -c ME_minn.f90 ${LIBS} 

main.o: main.f90 
	${FC} -c main.f90 ${LIBS} 

Hartree_Fock.o: Hartree_Fock.f90 
	${FC} -c Hartree_Fock.f90 ${LIBS} 

clean: 
	rm *.o
	rm *.mod
	rm ${p1} 
