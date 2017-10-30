#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
  int    i, n, n1 ;
  double p, p1, p2, pa, dp, sl,sg, kl, kg, p_ca, m ;
  FILE *fic ;
  fic = fopen("beton","w") ;
  p1 = 0. ;
  p2 = 2.e8 ;
  n  = 401 ;
  n1 = n-1 ;
  dp = (p2-p1)/n1 ;
  
  p_ca = 37547939.5 ;
  m    = 2.168403 ;
  for(i=0;i<n;i++) {
    p = p1 + i*dp ;
    sl = pow(1.+pow(p/p_ca,m/(m-1)),-1/m) ;
    kl = sqrt(sl)*pow(1. - pow(1. - pow(sl,m),1/m),2) ;
    sg = 1. - sl ;
    kg = sqrt(sg)*pow(1. - pow(sl,m),2/m) ;
    fprintf(fic,"%e %e %e %e\n",p,sl,kl,kg) ;
  }
  close(fic) ;
}
