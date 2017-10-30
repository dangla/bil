#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
int    i, n, n1 ;
double p, p1, p2, pa, dp, sl, kl ;
FILE *fic ;
fic = fopen("reze.tab1","w") ;
p1 = 0. ;
p2 = 1.e5 ;
n  = 401 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/0.5e4,2.2),-0.091) ;
  kl = pow(sl,16.2) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
fic = fopen("reze.tab2","w") ;
p1 = 0. ;
p2 = 1.e5 ;
n  = 401 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/0.13e4,2.45),-0.183) ;
  kl = pow(sl,6.2) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
fic = fopen("reze.tab3","w") ;
p1 = 0. ;
p2 = 1.e5 ;
n  = 401 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/0.38e4,2.26),-0.115) ;
  kl = pow(sl,11.2) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
}