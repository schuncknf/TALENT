# A Makefile with auto dep generation
# Note that the build will fail if header
# files are renamed. This is because the
# *.d dependency files will be out of date
# and point to non-existing header files
# We can solve this automatically but to
# keep the code nice and simple just do
# $> make clean (or rm -r dep/)
# $> make
# If header files have been renamed
#
# @author Camille Colle

CC=g++
WFLAGS=-Wall -Wextra -Wshadow #-Weffc++ # You can expect to see a buttload of warnings when enabling all these flags
CFLAGS= -O2 -std=c++11 -fPIC $(WFLAGS) `gsl-config --cflags`
LIBS=`gsl-config --libs` -lblas -larmadillo

TARGETS=quad_generator main test # files with a main()
OBJS=quadrature.o quad_ext.o# object files for building, no main()'s here !!!
SOLIBS=# shared object files go here
DEPDIR:=dep
DEPS=$(wildcard $(DEPDIR)/*.d)

.DEFAULT: all
all: $(TARGETS) $(SOLIBS)


-include $(DEPS)

$(SOLIBS): %.so: %.o
	@tput setaf 4; tput bold
	@printf "Building shared object $@ : $^ \n"
	@tput sgr0
	$(CC) -shared $< $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"

# dependency of each target is all the object files and it's own %.o file.
# Technically too much dependencies... (targets will generally not depend on all $(OBJS) files)
# But it's better to have too much dependencies then too few.
$(TARGETS): %: %.o $(OBJS)
	@tput setaf 1; tput bold
	@printf "Linking executable $@  :  $^ \n"
	@tput sgr0
	$(CC) $^ $(CFLAGS) $(INC) $(LIBS) -o $@
	@printf "\n"
	mv $@ ../run

%.o : %.cpp
	@echo DEPS INCLUDED ARE $(DEPS)
	@mkdir -p $(DEPDIR)
	@tput setf 6
	@printf "Building dependencies for $<\n"
	@tput sgr0
	$(CC) -MM $< $(CFLAGS) $(INC) > $(patsubst %.o,$(DEPDIR)/%.d,$@)
	@tput setf 6
	@printf "Building  $@ :  $^ \n"
	@tput sgr0
	$(CC) -c $< $(CFLAGS) $(INC) -o $@

.PHONY: clean

clean:
	rm -f *.o $(DEPDIR)/*.d *.so 
