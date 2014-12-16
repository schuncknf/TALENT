#ifndef SOURCE1_H_
#define SOURCE1_H_
#include <string>
using namespace std;
extern int Debug; 
extern string in_state_file;
extern string in_potential_file;
extern double (*Potential)(double ,double ,void* );

#define ST_n 0
#define ST_l 1
#define ST_m 2
#define ST_s 3
#define ST_t 4
#define ST_j 5

#define hbarc    197.327053
#define mass_p   938.27231
#define mass_n   939.56563
#define one_on_alpha  137.035989

#endif