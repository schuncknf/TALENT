#CC	= g++
CC	= gcc
CFLAGS	= -Wall
#LDFLAGS	= -L/usr/local/lib/ -lgsl -lgslcblas -lm
LDFLAGS	= -lgsl -lgslcblas -lm -Wno-deprecated 
EXEC	= run
DEPS	= param.h eigen.h solver.h potential.h gaulag.h sho.h check.h response.h
SRC	= main.c eigen.c solver.c potential.c gaulag.c sho.c check.c response.c

all:	$(EXEC)

run:	$(DEPS) $(SRC)
	$(CC) -o $@ $(SRC) $(CFLAGS) $(LDFLAGS)

clean:
	rm -rf *~

mrproper: clean
	rm -rf $(EXEC)

