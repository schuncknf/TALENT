

libquadrature.so: quad_ext.cpp quad_ext.hpp quadrature.cpp quadrature.hpp
	g++ -Wall -Wextra -Wshadow -O3 -fPIC -shared $^ -o $@

.PHONY: clean
clean:
	rm quadrature.so
