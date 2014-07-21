CC = icc
FC = gfortran
CFLAGS = -Wall -O3 -openmp -Wno-sign-compare -lgsl -lgslcblas
LDFLAGS = -llapack -lblas
MYINCLUDES = 
#LDFLAGS = -L /home/titus/Desktop/openblas -lopenblas
#MYINCLUDES = -I /home/titus/Desktop/openblas
OBJECTS = hydrogen.o functions.o harmonic_oscillator.o gauss_legendre.o





hydrogen : ${OBJECTS}
	${CC} ${CFLAGS} ${OBJECTS} ${MYINCLUDES} -o hydrogen ${LDFLAGS} -lgfortran

functions.o : functions.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c functions.cpp ${LDFLAGS}

hydrogen.o: hydrogen.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c hydrogen.cpp ${LDFLAGS}


harmonic_oscillator.o : harmonic_oscillator.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c harmonic_oscillator.cpp ${LDFLAGS}

gauss_legendre.o : gauss_legendre.c
	${CC} ${CFLAGS} ${MYINCLUDES} -c gauss_legendre.c ${LDFLAGS}









