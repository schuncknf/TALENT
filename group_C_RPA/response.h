#ifndef _pair_h
#define _pair_h

//#include "potential.h"

int npair;


typedef struct {
   int jpart;
   int lpart;
   int jlpart;
   int indpart;
   double epart;
   
   int jhole;
   int lhole;
   int jlhole;
   int indhole;
   double ehole;
   
   double ediff;
   
} pairs_t;

typedef struct{
   double energy;
   double *x;
   double *y;
} rpa_t;

//void create_pairs(waves_t*);

#endif
