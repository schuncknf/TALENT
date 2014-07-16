TARGET=hoFunction
OFILES=sphericalhofunc.o integratorGaussLegendre.o

CPPFLAGS=-Wall -g -std=c++0x -O3
MPICC=g++

INCLUDES=

LIBS= -larmadillo -lgsl



all: $(TARGET)


$(TARGET):  $(OFILES) $(TARGET).cpp
	$(MPICC) $(OFILES) $(TARGET).cpp -o $(TARGET) $(LIBS)

%.o: %.cpp %.h
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $< -c

clean:
	rm -f $(TARGET) $(OFILES)
