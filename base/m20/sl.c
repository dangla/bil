#include <stdio.h>
#include<math.h>

int main(void)
{
  double p,s,k,Pmin,Pmax,S1,S0,alpha,a,c;
  int    i,n;
  FILE   *fichier;
  
  
  Pmin = 500. ;
  Pmax = 1000. ;
  n = 2000 ;
  S0 = 0.09 ;
  S1 = 1. ;

  alpha = 3. ;
  
  
  a=1/(pow(S1,alpha)-pow(S0,alpha));
  c=-a*pow(S0,alpha);
  
  
  fichier=fopen("billes","w");
  
  for (p=Pmin;p<Pmax;p=p+((Pmax-Pmin)/(n-1))) {
    s=0.09+pow(pow(1-0.09,-57.53/56.53)+pow(p/784,57.53),(1/57.53)-1);
    
    if (s>S0)  k= a*pow(s,alpha)+c ;
    else k = 0 ;
    if (k>1) k=1;                   
    
    fprintf(fichier,"%e %e %e\n",p,s,k);
  }
  
  s=0.09+pow(pow(1-0.09,-57.53/56.53)+pow(Pmax/784,57.53),(1/57.53)-1);
  
  if (s>S0) 
    k=a*pow(s,alpha)+c;
  else k=0;
  if (k>1)
    k=1;  
  
  fprintf(fichier,"%e %e %e\n",Pmax,s,k);
  
  fclose(fichier);
  
  
  fichier=fopen("billes.s","w");
  
  for(i=0;i<n;i++) {
    s = S0 + i*(S1 - S0)/(n-1) ;
    if(s>S0) p = 784.*pow(pow(s-S0,-57.53/56.53) - pow(1.-S0,-57.53/56.53),1./57.53) ;
    else p = Pmax ;
    
    if (s>S0)  k=a*pow(s,alpha)+c;
    else k=0;
    if (k>1) k=1;                   
    
    fprintf(fichier,"%e %e %e\n",s,p,k);
  }
  
  fclose(fichier);
}
 
