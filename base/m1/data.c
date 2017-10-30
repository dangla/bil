#include <stdio.h>
#include<math.h>

int main(void)
 {
  double p,s,k,Pmin,Pmax,S1,S0,alpha,a,c;
  int  n;
  FILE *fichier;
  
  
  printf("valeur de Pmin:");
  scanf("%lf", &Pmin);
  printf("valeur de Pmax:");
  scanf("%lf", &Pmax);
  printf("nombre de points :");
  scanf("%d", &n);
  printf("coeff alpha, de Kr(S) :");
  scanf("%lf", &alpha);
  printf("S0 de Kr(S) :");
  scanf("%lf", &S0);
  printf("S1 de Kr(S) :");
  scanf("%lf", &S1);
  
  
  a=1/(pow(S1,alpha)-pow(S0,alpha));
  c=-a*pow(S0,alpha);
  
  
  fichier=fopen("billes","w");
    
  for (p=Pmin;p<Pmax;p=p+((Pmax-Pmin)/(n-1)))
  {
    s=0.09+pow(pow(1-0.09,-57.53/56.53)+pow(p/784,57.53),(1/57.53)-1);
    
    if (s>S0) 
       k=a*pow(s,alpha)+c;
       else k=0;
    if (k>1)
       k=1;                   
          
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
 }
 
