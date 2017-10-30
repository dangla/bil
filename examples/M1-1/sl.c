#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
int    i, n, n1 ;
double p, p1, p2, pa, dp, sl, kl, kg ;
FILE *fic ;
fic = fopen("tab1","w") ;
p1 = 0. ;
p2 = 5.e8 ;
n  = 401 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/1.5e6,1.06383),-0.06) ;
  kl = pow(1.+pow(p/3.e6,2.),-0.5) ;
  kg = (1.-sl)*(1.-sl)*(1.-pow(sl,5./3.)) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
fic = fopen("tab2","w") ;
p1 = 0. ;
p2 = 1.e7 ;
n  = 201 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/10.e6,1.7),-0.4117) ;
  kl = pow(1.+pow(p/10.e6,2.),-1.) ;
  kg = (1.-sl)*(1.-sl)*(1.-pow(sl,5./3.)) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
}