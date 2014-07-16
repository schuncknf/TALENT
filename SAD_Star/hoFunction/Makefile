#TARGET=hoFunction vmatrixgenerator.o
TARGET=vmatrixGeneratorC
OFILES=sphericalhofunc.o integratorGaussLegendre.o 

CPPFLAGS=-Wall -g -std=c++0x -O3
MPICC=g++

INCLUDES=

LIBS= -g -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -larmadillo


all: $(TARGET)


$(TARGET):  $(OFILES) $(TARGET).cpp
	$(MPICC) $(OFILES) $(TARGET).cpp -o $(TARGET) $(LIBS)

%.o: %.cpp %.h
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $< -c
clean:
	rm -f $(TARGET) $(OFILES)
