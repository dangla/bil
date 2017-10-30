#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
int    i, n, n1 ;
 double p, p1, p2, pa, dp, sl, kl ;
FILE *fic ;
fic = fopen("roche","w") ;
p1 = 0. ;
p2 = 2.e8 ;
n  = 401 ;
n1 = n-1 ;
dp = (p2-p1)/n1 ;
for(i=0;i<n;i++)
  {
  p = p1 + i*dp ;
  sl = pow(1.+pow(p/1.e7,1./(1.-0.412)),-0.412) ;
  kl = 1./(1.+pow(pow(sl,-2.429)-1.,1.176)) ;
  fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
close(fic) ;
}
