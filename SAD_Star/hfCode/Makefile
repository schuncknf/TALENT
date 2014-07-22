TARGET=main
TESTTARGET= testMain
OFILES=sphericalhofunc.o integratorGaussLegendre.o integratorGaussLaguerre.o integratorGaussLegendreGSL.o VCoulombMatrixGenerator.o

CPPFLAGS=-Wall -g -std=c++0x -O3
MPICC=g++

INCLUDES= -I/opt/local/include

LIBS= -lgsl -lgslcblas -L/opt/local/lib -larmadillo
TESTLIB= $(LIBS) -lboost_unit_test_framework


all: $(TARGET)

$(TARGET):  $(OFILES) $(TARGET).cpp
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $(OFILES) $(TARGET).cpp -o $(TARGET) $(LIBS)

%.o: %.cpp %.h
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $< -c

tests:  $(OFILES) $(TESTTARGET).cpp
	$(MPICC) $(CPPFLAGS) $(INCLUDES) $(OFILES) $(TESTTARGET).cpp -o $(TESTTARGET) $(TESTLIB)

clean:
	rm -f $(TARGET) $(TESTTARGET) $(OFILES)
