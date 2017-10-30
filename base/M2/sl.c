#include <stdio.h>
#include <math.h>
int main()
{
  int    i,n ;
  double sl,p,p1,p2,dp,kl,kg,sl_e ;
  FILE   *fic ;

  fic = fopen("sol","w") ;
  p1 = 0. ;
  p2 = 1.5e4 ;
  n  = 101 ;
  dp = (p2-p1)/(n-1) ;
  for(i=0;i<n;i++) {
      p = p1 + i*dp ;
      sl = 1. - 0.096863*pow(p/9.80665e3,2.4279495) ;
      kl = 1. - 2.207*pow(1. - sl,1.0121) ;
      sl_e = (sl - 0.2)/0.8 ;
      if(sl_e > 0.) kg = (1. - sl_e)*(1. - sl_e)*(1. - pow(sl_e,5./3.)) ; else kg = 1. ;
      fprintf(fic,"%e %e %e %e\n",p,sl,kl,kg) ;
    }
}
