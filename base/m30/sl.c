#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(void)
{
  int    i, n = 101, n1 = n - 1 ;
  double p, p1 = 1.e-1, p2 = 1.e6, pa, dp = (p2 - p1)/n1, sl, kl, se, sr, nl, ml ;
  FILE *fic ;
  fic = fopen("fracture_log","w") ;
  sr = 0. ;
  nl = 2. ;
  ml = 1. - 1./nl ;
  for(i=0;i<n;i++)
    {
      p = p1 + i*dp ;
      p = p1*pow(p2/p1,((double) i)/n1) ;
      sl = sr + (1. - sr)*pow(1.+pow(p/9.81e2,nl),-ml) ;
      se = (sl - sr)/(1 - sr) ;
      kl = pow(se,0.5)*pow(1. - pow(1. - pow(se,1./ml),ml),2.) ;
      fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
  close(fic) ;
  fic = fopen("matrice_log","w") ;
  sr = 0.2105 ;
  nl = 1.5 ;
  ml = 1. - 1./nl ;
  for(i=0;i<n;i++)
    {
      p = p1 + i*dp ;
      p = p1*pow(p2/p1,((double) i)/n1) ;
      sl = sr + (1. - sr)*pow(1.+pow(p/1.96e4,nl),-ml) ;
      se = (sl - sr)/(1 - sr) ;
      kl = pow(se,0.5)*pow(1. - pow(1. - pow(se,1./ml),ml),2.) ;
      fprintf(fic,"%e %e %e\n",p,sl,kl) ;
  }
  close(fic) ;
  fic = fopen("interface_log","w") ;
  sr = 0.2105 ;
  nl = 1.5 ;
  ml = 1. - 1./nl ;
  for(i=0;i<n;i++)
    {
      p = p1 + i*dp ;
      p = p1*pow(p2/p1,((double) i)/n1) ;
      sl = sr + (1. - sr)*pow(1.+pow(p/1.96e4,nl),-ml) ;
      se = (sl - sr)/(1 - sr) ;
      kl = pow(se,0.5)*pow(1. - pow(1. - pow(se,1./ml),ml),2.) ;
      fprintf(fic,"%e %e\n",p,kl) ;
    }
  close(fic) ;
}
