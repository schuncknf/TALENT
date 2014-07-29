

all: quad_generator libquadrature.so

quad_generator: quad_generator.cpp quadrature.hpp quadrature.cpp quad_ext.hpp quad_ext.cpp
	g++ -Wall -Wextra -Wshadow -O3 -fPIC $^ -o $@ 

libquadrature.so: quad_ext.cpp quad_ext.hpp quadrature.cpp quadrature.hpp
	g++ -Wall -Wextra -Wshadow -O3 -fPIC -shared $^ -o $@

.PHONY: clean
clean:
	rm quadrature.so
