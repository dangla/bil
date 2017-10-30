#include <stdio.h>
#include <math.h>

/* doit etre compiler avec m44.c ou la version de m150.c de TQN */

extern double activite(double,double,double,double,double,int,double*) ;

int main(void)
{
  int    n = 100,n1 = n - 1 ;
  double M_na = 22.95e-3,M_cl = 35.45e-3,M_nacl = M_na + M_cl ;
  double c1 = 0,c2 = 350/M_nacl,dc = (c2 - c1)/n1 ;
  double c_cl,c_oh = 1.e-10,c_na,c_k = c_oh,T = 293. ;
  double lna[4],lna_w ;
  int    i ;

  for(i=0;i<n;i++) {
    c_cl = c1 + i*dc ;
    c_na = c_cl ;
    c_k = c_oh ; /* electroneutralite */
    lna_w = activite(c_cl,c_oh,c_na,c_k,T,1,lna) ;
    printf("%e %e %e %e\n",c_cl*M_nacl,lna_w,lna[0],lna[2]) ;
  }
}
