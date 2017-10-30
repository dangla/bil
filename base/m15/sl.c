#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../../src/DataSet/Functions.h"

int main(void)
{
  int    i, n, n1 ;
  double p, p1, p2, pa, dp, sl, kl, l ;
  double x[6] = {0,10.e3,20.e3,40.e3,80.e3,160.e3} ;
  double f[6] = {1,1,2.22222,3,4.0555555,6.666666} ;
  Function_t fn = {6,x,f} ;
  FILE *fic ;
  fic = fopen("sterrebeek","w") ;
  p1 = 10.e3 ;
  p2 = 210.e3 ;
  n  = 101 ;
  n1 = n-1 ;
  dp = (p2-p1)/n1 ;
  for(i=0;i<n;i++)
    {
      p = p1 + i*dp ;
      sl = pow(p/10.e3,-1./2.5) ;
      kl = pow(1.+pow(p/3.e6,2.),-0.5) ;
      l  = Function_ComputeValue(&fn,p) ;
      fprintf(fic,"%e %e %e %e\n",p,sl,kl,l) ;
    }
  close(fic) ;
}
