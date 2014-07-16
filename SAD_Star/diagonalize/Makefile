TARGET=diagonalize

CPPFLAGS=-Wall -g -std=c++0x -O3
MPICC=g++



INCLUDES=


LIBS= -larmadillo
#-llapack


all: 
	$(MPICC) $(INCLUDES) $(TARGET).cpp -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET) 
