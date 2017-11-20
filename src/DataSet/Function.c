#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Function.h"




Function_t*  (Function_New)(const int n)
{
  Function_t* function = (Function_t*) malloc(sizeof(Function_t)) ;
    
  if(!function) {
    arret("Function_New(1)") ;
  }
    
  Function_GetNbOfPoints(function) = n ;

  {
    double* tm = (double*) malloc(2*n*sizeof(double)) ;
    
    if(!tm) {
      arret("Function_New(2)") ;
    }
    
    Function_GetXValue(function) = tm ;
    Function_GetFValue(function) = tm + n ;
  }
  
  return(function) ;
}



double (Function_ComputeValue)(Function_t* fn,double t)
{
  if(fn) {
    int   nb_points = Function_GetNbOfPoints(fn) ;
    double* tm = Function_GetXValue(fn) ;
    double* ft = Function_GetFValue(fn) ;
    
    /*
      Cas t < t[0] 
    */
    if(t <= tm[0]) {
      return(ft[0]) ;
      
    /*
      Cas t > t[n-1] 
    */
    } else if(t >= tm[nb_points - 1]) {
      return(ft[nb_points - 1]) ;
      
    /*
      Cas t[i1] <= t <= t[i2]
      Calcul des deux points i1=Min(i) et i2=Max(i) du tableau fn
      correspondant au plus petit intervalle de temps [t1;t2]
      contenant t. Deux cas de figure se presentent:
      1) t2 > t1 et i2 = i1+1;
      2) t2 = t1 et i2 >= i1 avec i2-i1 maximum.
    */
    } else {
      int    i1 = 0, i2 = nb_points - 1 ;
      double t1, t2, f1, f2 ;
      int    i ;
      
      for(i = 0 ; i < nb_points ; i++) {
        if(t >= tm[i]) i1 = i ;
        if(t <= tm[i]) break  ;
      }
      
      for(i = i1 ; i < nb_points ; i++) if(t < tm[i]) break ;
      
      for( ; i >= 0 ; i--) {
        if(t <= tm[i]) i2 = i ;
        if(t >= tm[i]) break  ;
      }
      
      t1 = tm[i1] ;
      f1 = ft[i1] ;
      t2 = tm[i2] ;
      f2 = ft[i2] ;
      
      /* 1) Intervalle non nul */
      if(t2 > t1) {
        return (f1 + (f2 - f1)*(t - t1)/(t2 - t1)) ;
        
      /* 2) Intervalle nul */
      } else {
        /* un seul point */
        if(i1 == i2) {
          return (f1) ;
        /* plusieurs points */
        } else {
          arret("plusieurs pas de temps nuls (fonction)") ;
          return(0.) ;
        }
      }
    }
  }
  
  return(0.) ;
}
