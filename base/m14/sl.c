#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
int    i, n, n1 ;
double p, p1, p2, pa, dp, sl, kl, kg ;
FILE *fic ;
fic = fopen("clay","w") ;
p1 = 0. ;
p2 = 1.e9 ;
n  = 1001 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/1.8e7,1.613),-0.38) ;
  kl = pow(1.+pow(p/3.e6,2.),-0.5) ;
  kg = (1.-sl)*(1.-sl)*(1.-pow(sl,5./3.)) ;
  fprintf(fic,"%e %e %e %e\n",p,sl,kl,kg) ;
  }
close(fic) ;
}
