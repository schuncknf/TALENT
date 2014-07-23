CC = icc
FC = gfortran
CFLAGS = -Wall -O3 -openmp -Wno-sign-compare -lgsl -lgslcblas
LDFLAGS = -llapack -lblas
MYINCLUDES = 
#LDFLAGS = -L /home/titus/Desktop/openblas -lopenblas
#MYINCLUDES = -I /home/titus/Desktop/openblas
OBJECTS = hf.o functions.o harmonic_oscillator.o gauss_legendre.o minnesotame.o hffunctions.o





hf : ${OBJECTS}
	${CC} ${CFLAGS} ${OBJECTS} ${MYINCLUDES} -o hf ${LDFLAGS} -lgfortran

functions.o : functions.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c functions.cpp ${LDFLAGS}

minnesotame.o : minnesotame.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c minnesotame.cpp ${LDFLAGS}

hffunctions.o : hffunctions.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c hffunctions.cpp ${LDFLAGS}

hf.o: hf.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c hf.cpp ${LDFLAGS}


harmonic_oscillator.o : harmonic_oscillator.cpp
	${CC} ${CFLAGS} ${MYINCLUDES} -c harmonic_oscillator.cpp ${LDFLAGS}

gauss_legendre.o : gauss_legendre.c
	${CC} ${CFLAGS} ${MYINCLUDES} -c gauss_legendre.c ${LDFLAGS}









