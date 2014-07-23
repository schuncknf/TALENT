CC	= gcc
CFLAGS	= -O3 -Wall -ansi -pedantic
LDFLAGS	= -L/usr/local/lib/ -lgsl -lgslcblas -lm
#LDFLAGS	= -llapacke -llapack -lblas -lm -Wno-deprecated 
EXEC	= run
DEPS	= param.h eigen.h solver.h potential.h gaulag.h sho.h
#DEPS	= param.h eigen.h solver.h potential.h gaulag.h lapacke.h sho.h
SRC	= main.cc eigen.cc solver.cc potential.cc gaulag.cc sho.cc
OBJ	= $(SRC:.c=o.)

all:	$(EXEC)

run:	$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

%.o: %.c $(DEPS)
	$(cc) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o *~

mrproper: clean
	rm -rf $(EXEC)

