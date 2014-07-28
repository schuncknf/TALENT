#include "quadrature.hpp"

/** @param kind_ the kind of quadrature, 5 is generalized Gauss Laguerre
 *  @param alpha_ alpha is zero for plain Gauss Laguerre
 *  @param beta_  beta does not matter for Gauss Laguerre
 *  @param nt_  number of knots to write
 *  for more info see the file quad_ext.hpp
 *  btw the underscores are just to prevent -Wshadow warnings
*/
Quad::Quad(const int kind_,const double alpha_, const double beta_, const int nt_) : kind(kind_),alpha(alpha_),beta(beta_),nt(nt_)
{
	knots = new double[nt];
	weights = new double[nt];
	quad::cdgqf(nt,kind,alpha,beta,knots,weights); /**< HERE is the call to the code i got from the internet, calculates knots and weights */
}

Quad::Quad(const Quad& rhs){ // copy constructor to make sure the dynamically allocated knots and weights are safely copied over
	kind  = rhs.kind;
	alpha = rhs.alpha;
	beta  = rhs.beta;
	nt    = rhs.nt;
	
	knots = new double[nt];
	weights = new double[nt];
	for (int i=0; i<nt; i++){
		knots[i]  = rhs.knots[i];
		weights[i]= rhs.weights[i];
	}
}

Quad& Quad::operator=(const Quad& rhs){
	if (knots)
		delete[] knots;
	if (weights)
		delete[] weights;

	kind  = rhs.kind;
	alpha = rhs.alpha;
	beta  = rhs.beta;
	nt    = rhs.nt;
	
	knots = new double[nt];
	weights = new double[nt];
	for (int i=0; i<nt; i++){
		knots[i]  = rhs.knots[i];
		weights[i]= rhs.weights[i];
	}
	return *this;
}

Quad::Quad() : kind(-1),alpha(-1),beta(-1),nt(-1),knots(NULL),weights(NULL) {}

Quad::Quad(const char* inputfile) : kind(-1), alpha(-1), beta(-1),nt(-1),knots(NULL),weights(NULL) {
	read(inputfile);
}

Quad::~Quad(){
	if (knots != NULL){
		delete[] knots;
		delete[] weights;
	}
}

void Quad::write(const char* outputfile){
	std::ofstream file(outputfile,std::ios::binary);
	if (file.is_open()){
		file.write((char*) &kind,sizeof(int));
		file.write((char*) &alpha,sizeof(double));
		file.write((char*) &beta,sizeof(double));
		file.write((char*) &nt,sizeof(int));
		for (int i=0; i<nt; i++){
			file.write((char*) &knots[i],sizeof(double));
			file.write((char*) &weights[i],sizeof(double));
		}
		file.close();
	} else {
		std::cerr << "[Error] Unable to open file : " << outputfile << std::endl;
		exit(-1);
	}
}

void Quad::read(const char* inputfile){
	std::ifstream file(inputfile,std::ios::binary);	
	if (file.is_open()){
		file.read((char*) &kind, sizeof(int));
		file.read((char*) &alpha,sizeof(double));
		file.read((char*) &beta, sizeof(double));
		file.read((char*) &nt, sizeof(int));
		weights = new double[nt];
		knots   = new double[nt];
		for (int i=0; i<nt; i++){
			file.read((char*) &knots[i],sizeof(double));
			file.read((char*) &weights[i],sizeof(double));
		}			
	} else {
		std::cerr << "[Error] Unable to open file : " << inputfile << std::endl;
		exit(-1);
	}
}

double Quad::integrate( double (*f)(double,void*), void* f_params){ // Gauss Laguerre Quadrature Integration
	assert(kind==5); // only Gauss Laguerre Quadrature supported for now, can be extended
	double res = 0;
	for (int i=0; i<nt; i++) {
		/** exp(knots[i]) can become numerically unstable if knots become large, instead of multiplying sum logs and exp later!
		* a*b = e^ (ln(a) + ln(b) )
		* original formula that was found to be numerically unstable given below
		* note that f[i] can be negative so some sign checking is done to expand f(t) as exp(ln(f(t))) if pos and as -exp(ln(|f(t)|)) if neg
		* res += weights[i]*f(knots[i],f_params)*exp(knots[i])*pow(knots[i],-alpha); // exp(knots) needed to compensate function see wiki for expl.
		*/
		const double ft = f(knots[i],f_params);
		if (fabs(ft) > 0.) {
			const int sign  = std::signbit( ft ) ? -1 : 1;
			res += sign*exp( log( weights[i]*fabs(ft) ) + knots[i] - alpha*log(knots[i]) );
		}
	}
	return res;
}

double Quad::integrate( double (*f)(double,double,void*), void* f_params){
	assert(kind==5);
	double res=0;
	for (int i=0; i<nt;i++){
		for (int j=0;j<nt;j++){
			const double ft = f(knots[i],knots[j],f_params);
			if (fabs(ft)>0.){
				const int sign = std::signbit(ft) ? -1 : 1;
				/**
				 * over/under flow unsafe version
				 * res+= weights[i]*weights[j]*f(knots[i],f_params)*exp(knots[i])*pow(knots[i],-alpha)*exp(knots[j])*pow(knots[j],-alpha)
				*/
				res += sign*exp( log( weights[i]*weights[j]*fabs(ft) ) + knots[i] - alpha*log(knots[i]) + knots[j] - alpha*log(knots[j]) );
			}
		}
	}
	return res;
}
