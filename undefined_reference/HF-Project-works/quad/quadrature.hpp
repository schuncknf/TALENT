/** Quadrature Integration
* for now only Gauss Laguerre supported although more to be found
* in the file quad_ext.hpp
*/

#ifndef QUADRATURES_HPP
#define QUADRATURES_HPP
#include "quad_ext.hpp"
#include <cassert>

class Quad{
	public:
		Quad();
		Quad(const int kind, const double alpha, const double beta, const int nt);		
		Quad(const char*); // construct quadrature from file
		Quad(const Quad&); // copy constructor
		Quad& operator= (const Quad&); // assignment operator
		~Quad();
		void write(const char* output);
		void read(const char* input);
		double integrate( double (*f)(double,void*), void* f_params); // 1D Quadrature Integration
		double integrate( double (*f)(double,double,void*), void* f_params); // 2D Quadrature Integration
	private:
		int kind;
		double alpha,beta;
		int nt;
		double *knots,*weights;
}; // class Quad
#endif // QUADRATURES_HPP
