#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
  int    i, n, n1 ;
  double p, p1, p2, pa, dp, sl, kl, p_ca, m, b, tau_l, tau_g ;
  FILE *fic ;
  fic = fopen("beton","w") ;
  p1 = 0. ;
  p2 = 1.e10 ;
  n  = 4001 ;
  n1 = n-1 ;
  dp = (p2-p1)/n1 ;
  b    = 3.2 ;
  
  p_ca = 1.86e7 ;
  m    = 2.27 ;
  for(i=0;i<n;i++) {
    p = p1 + i*dp ;
    sl = pow(1.+pow(p/p_ca,m/(m-1)),-1/m) ;
    kl = sqrt(sl)*pow(1. - pow(1. - pow(sl,m),1/m),2) ;
    tau_l = 1./(sl*(1+625*pow(1-sl,4))) ;
    tau_g = pow(1-sl,b) ;
    fprintf(fic,"%e %e %e %e %e\n",p,sl,kl,tau_l,tau_g) ;
  }
  close(fic) ;
}
