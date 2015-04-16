#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <R.h>
 
void gibbsGSL(int *Np,int *thinp,int *seedp,double *xvec,double *yvec)
{
  int i,j;
  int N=*Np,thin=*thinp,seed=*seedp;
  gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);
  double x=0;
  double y=0;
  for (i=0;i<N;i++) {
    for (j=0;j<thin;j++) {
      x=gsl_ran_gamma(r,3.0,1.0/(y*y+4));
      y=1.0/(x+1)+gsl_ran_gaussian(r,1.0/sqrt(2*x+2));
    }
    xvec[i]=x; yvec[i]=y;
  }
}