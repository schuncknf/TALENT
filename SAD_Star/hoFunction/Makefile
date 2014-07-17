#TARGET=hoFunction vmatrixgenerator.o
TARGET=vmatrixGeneratorC
OFILES=sphericalhofunc.o integratorGaussLegendre.o 

CPPFLAGS=-Wall -g -std=c++0x -O3
MPICC=g++

INCLUDES= -I/opt/local/include

LIBS= -lgsl -lgslcblas -L/opt/local/lib -larmadillo


all: $(TARGET)


$(TARGET):  $(OFILES) $(TARGET).cpp
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $(OFILES) $(TARGET).cpp -o $(TARGET) $(LIBS)

tests:  $(OFILES) tests.cpp
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $(OFILES) tests.cpp -o tests $(LIBS)

%.o: %.cpp %.h
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $< -c
clean:
	rm -f $(TARGET) $(OFILES)
