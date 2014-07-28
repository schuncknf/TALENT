#CC	= g++
CC	= gcc
CFLAGS	= -Wall -O2
#LDFLAGS	= -L/usr/local/lib/ -lgsl -lgslcblas -lm
LDFLAGS	= -lgsl -lgslcblas -lm -Wno-deprecated 
EXEC	= run
DEPS	= param.h eigen.h solver.h potential.h gaulag.h sho.h
SRC	= main.c eigen.c solver.c potential.c gaulag.c sho.c

all:	$(EXEC)

run:	$(DEPS) $(SRC)
	$(CC) -o $@ $(SRC) $(CFLAGS) $(LDFLAGS)

clean:
	rm -rf *~

mrproper: clean
	rm -rf $(EXEC)

