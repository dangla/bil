#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <strings.h>
#include <stdarg.h>

#include "Message.h"
#include "Mry.h"
#include "Tools/Math.h"
#include "Curve.h"


static double courbe_nor(double,Curve_t*) ;
static double dcourbe_nor(double,Curve_t*) ;
static double icourbe_nor(double,Curve_t*) ;
static double courbe_log(double,Curve_t*) ;
static double dcourbe_log(double,Curve_t*) ;
static double icourbe_log(double,Curve_t*) ;


/* Extern functions */

Curve_t* (Curve_Create)(unsigned int n_points)
{
  Curve_t* curve   = (Curve_t*) Mry_New(Curve_t) ;

  Curve_GetNbOfPoints(curve) = n_points ;
  
  
  /* Allocate memory space for the values */
  {
    double* mry = (double *) Mry_New(double[n_points + 2]) ;
    
    Curve_GetXRange(curve) = mry ;
    Curve_GetYValue(curve) = mry + 2 ;
  }
  
  
  /* Allocate memory space for the names of axis */
  {
    char* name = (char*) Mry_New(char[2*Curve_MaxLengthOfCurveName]) ;
    
    Curve_GetNameOfXAxis(curve) = name ;
    Curve_GetNameOfYAxis(curve) = name + Curve_MaxLengthOfCurveName ;
    
    name[0] = '\0' ;
    name[Curve_MaxLengthOfCurveName] = '\0' ;
  }
  
  return(curve) ;
}



void (Curve_Delete)(void* self)
{
  Curve_t* curve   = (Curve_t*) self ;
  
  {
    double* mry = Curve_GetXRange(curve) ;
    
    if(mry) {
      free(mry) ;
      Curve_GetXRange(curve) = NULL ;
    }
  }
  
  {
    char* name = Curve_GetNameOfXAxis(curve) ;
    
    if(name) {
      free(name) ;
      Curve_GetNameOfXAxis(curve) = NULL ;
    }
  }
}



Curve_t* Curve_CreateDerivative(Curve_t* curve)
{
  int n_points    = Curve_GetNbOfPoints(curve) ;
  Curve_t* dcurve = Curve_Create(n_points) ;
  double* dy      = Curve_GetYValue(dcurve) ;
  double* x       = Curve_CreateSamplingOfX(curve) ;
  int i ;
  
  Curve_GetNbOfPoints(dcurve) = Curve_GetNbOfPoints(curve) ;
  Curve_GetXRange(dcurve)[0]  = Curve_GetXRange(curve)[0] ;
  Curve_GetXRange(dcurve)[1]  = Curve_GetXRange(curve)[1] ;
  Curve_GetScaleType(dcurve)  = Curve_GetScaleType(curve) ;
  
  for(i = 0 ; i < n_points ; i++) {
    dy[i] = Curve_ComputeDerivative(curve,x[i]) ;
  }
  
  free(x) ;

  return(dcurve) ;
}



Curve_t* Curve_CreateIntegral(Curve_t* curve)
{
  int n_points    = Curve_GetNbOfPoints(curve) ;
  Curve_t* icurve = Curve_Create(n_points) ;
  double* intydx  = Curve_GetYValue(icurve) ;
  double* x       = Curve_CreateSamplingOfX(curve) ;
  int i ;
  
  Curve_GetNbOfPoints(icurve) = Curve_GetNbOfPoints(curve) ;
  Curve_GetXRange(icurve)[0]  = Curve_GetXRange(curve)[0] ;
  Curve_GetXRange(icurve)[1]  = Curve_GetXRange(curve)[1] ;
  Curve_GetScaleType(icurve)  = Curve_GetScaleType(curve) ;
  
  for(i = 0 ; i < n_points ; i++) {
    intydx[i] = Curve_ComputeIntegral(curve,x[i]) ;
  }
  
  free(x) ;

  return(icurve) ;
}



Curve_t* Curve_CreateInverse(Curve_t* curve,const char scale)
/** Create the inverse of a function assumed as monotoneous.
 *  The sampling is defined by scale.
 *  Return a pointer to Curve_t.
 */
{
  int     n = Curve_GetNbOfPoints(curve) ;
  Curve_t* evruc = Curve_Create(n) ;
  
  Curve_GetNbOfPoints(evruc) = Curve_GetNbOfPoints(curve) ;
  Curve_GetScaleType(evruc)  = scale ;
  Curve_GetXRange(evruc)[0]  = Curve_GetYValue(curve)[0] ;
  Curve_GetXRange(evruc)[1]  = Curve_GetYValue(curve)[n - 1] ;
  
  {
    double *xo  = Curve_CreateSamplingOfX(curve) ;
    double *yo  = Curve_GetYValue(curve) ;
    double *yi  = Curve_CreateSamplingOfX(evruc) ;
    double *xi  = Curve_GetYValue(evruc) ;
    double ymax = MAX(fabs(yo[0]),fabs(yo[n - 1])) ;
    double err ;
    double tol = 1.e-10*ymax ;
    int i ;
    
    /* Compute the xi's */
    xi[0    ] = Curve_GetXRange(curve)[0] ;
    xi[n - 1] = Curve_GetXRange(curve)[1] ;
    
    for(i = 1 ; i < n - 1 ; i++) {
      double y = yi[i] ;
      double x = xi[i - 1] ;
      
      /* Initialize x more efficiently */
      {
        int j = 1 ;
        
        while( (y - yo[j - 1])*(y - yo[j]) > 0) j++ ;
        
        x = xo[j - 1] ;
      }
      
      /* Given y, find x = Inv[f](y) i.e. y = f(x) */
      {
        int    j = 0 ;
    
        do {
          double f  = Curve_ComputeValue(curve,x) ;
          double df = Curve_ComputeDerivative(curve,x) ;
          double dx ;
      
          if(df == 0) {
            arret("Curve_CreateInverse: inversion impossible") ;
          }
          
          dx = - (f - y)/df ;
      
          x += dx ;
      
          err = fabs(f - y) ;
      
          if(j++ > 50) {
            printf("f    = %e\n",f) ;
            printf("df   = %e\n",df) ;
            printf("dx   = %e\n",dx) ;
            arret("Curve_CreateInverse: no convergence") ;
          }
        } while(err > tol) ;
        
        xi[i] = x ;
      }
    }
    
    free(yi) ;
    free(xo) ;
  }
  
  return(evruc) ;
}



double* Curve_CreateSamplingOfX(Curve_t* curve)
{
  int n_points = Curve_GetNbOfPoints(curve) ;
  double* x    = Curve_GetXRange(curve) ;
  char scale   = Curve_GetScaleType(curve) ;
  int ni    = n_points - 1 ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double* xs = (double*) Mry_New(double[n_points]) ;
  
  if(scale == 'n') {
    double dx = (x2 - x1) ;
    int i ;
    
    for(i = 0 ; i < n_points ; i++) {
      double n = ((double) i)/ni ;
      
      xs[i] = x1 + n*dx ;
    }
    
  } else if(scale == 'l') {
    double ratio = x2/x1 ;
    int i ;
    
    if(x1 <= 0. || x2 <= 0.) {
      arret("Curve_CreateSamplingOfX(1)") ;
    }
    
    for(i = 0 ; i < n_points ; i++) {
      double n = ((double) i)/ni ;
      
      xs[i] = x1*pow(ratio,n) ;
    }
    
  } else {
    arret("Curve_CreateSamplingOfX(2)") ;
  }

  return(xs) ;
}



double Curve_ComputeValue(Curve_t* cb,double a)
/** Return the value at a */
{
  if(cb) {
    if(Curve_GetScaleType(cb) == 'n') return(courbe_nor(a,cb)) ;
    else if(Curve_GetScaleType(cb) == 'l') return(courbe_log(a,cb)) ;
    else arret("Curve_ComputeValue: option non prevue") ;
  } else {
    arret("Curve_ComputeValue: undefined curve") ;
  }
  
  return(0.) ;
}



double Curve_ComputeDerivative(Curve_t* cb,double a)
/** Return the derivative at a */
{
  if(cb) {
    if(Curve_GetScaleType(cb) == 'n') return(dcourbe_nor(a,cb)) ;
    else if(Curve_GetScaleType(cb) == 'l') return(dcourbe_log(a,cb)) ;
    else arret("Curve_ComputeDerivative: option non prevue") ;
  } else {
    arret("Curve_ComputeDerivative: undefined curve") ;
  }
  
  return(0.) ;
}



double Curve_ComputeIntegral(Curve_t* cb,double a)
/** Return the integral from begin to a */
{
  if(cb) {
    if(Curve_GetScaleType(cb) == 'n') return(icourbe_nor(a,cb)) ;
    else if(Curve_GetScaleType(cb) == 'l') return(icourbe_log(a,cb)) ;
    else arret("Curve_ComputeIntegral: option non prevue") ;
  } else {
    arret("Curve_ComputeIntegral: undefined curve") ;
  }
  
  return(0.) ;
}



char* Curve_PrintInFile(Curve_t* curve)
{
  int     n = Curve_GetNbOfPoints(curve) ;
  double* x = Curve_CreateSamplingOfX(curve) ;
  double* y = Curve_GetYValue(curve) ;
  /* Create a temporary filename */
  /* char*   targetfile = (char*) tempnam(".",NULL) ; *//* no portability */
  char*   targetfile = (char*) tmpnam(NULL) ;
  FILE*   target     = fopen(targetfile,"w") ;
  int i ;
  
  for(i = 0 ; i < n ; i++) {
    fprintf(target,"%e %e\n",x[i],y[i]) ;
  }
  
  fclose(target) ;
  
  free(x) ;
  
  return(targetfile) ;
}




/* Intern functions */

double courbe_nor(double a,Curve_t* cb)
{
  /* Retourne la valeur de la courbe en a */
  int    ni = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0] ;
  double a2 = Curve_GetXRange(cb)[1] ;

  if(a <= a1) return(Curve_GetYValue(cb)[0]) ;
  else if(a >= a2) return(Curve_GetYValue(cb)[ni]) ;
  else {
    double da = (a2 - a1)/ni ;
    double r  = (a - a1)/da ;
    int i  = floor(r) ;
    double a0 = a1 + i*da ;
    double dv = (Curve_GetYValue(cb)[i+1] - Curve_GetYValue(cb)[i])/da ;
    double v  = Curve_GetYValue(cb)[i] + dv*(a - a0) ;
    return(v) ;
  }
}



double dcourbe_nor(double a,Curve_t* cb)
{
  /* Retourne la derivee de la courbe en a */
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double da  = (Curve_GetXRange(cb)[1] - Curve_GetXRange(cb)[0])/n_i ;

  return((courbe_nor(a + da,cb) - courbe_nor(a - da,cb))*0.5/da) ;
}



double icourbe_nor(double a,Curve_t* cb)
/* Return the integral computed from cb */ 
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double *x = Curve_GetXRange(cb) ;
  double *y = Curve_GetYValue(cb) ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double dx = (x2 - x1)/n_i ;
  double fa = Curve_ComputeValue(cb,a) ;
  double intydx = 0 ;
  int    i ;
  
  for(i = 0 ; x1 + (i + 1)*dx < a ; i++) {
    intydx += y[i] + y[i + 1] ;
  }
  
  intydx *= dx*0.5 ;
  
  intydx += (y[i] + fa)*(a - x1 - i*dx)*0.5 ;
  
  return(intydx) ;
}



double courbe_log(double a,Curve_t* cb)
/* Retourne la valeur en a de la courbe echantillonnee en base log10 */
{
  int    ni = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0] ;
  double a2 = Curve_GetXRange(cb)[1] ;

  if(a <= a1) return(Curve_GetYValue(cb)[0]) ;
  else if(a >= a2) return(Curve_GetYValue(cb)[ni]) ;
  else {
    double loga1 = log10(a1) ;
    double loga2 = log10(a2) ;
    double dloga = (loga2 - loga1)/ni ;
    double loga  = log10(a) ;
    double r  = (loga - loga1)/dloga ;
    int    i  = floor(r) ;
    if(i >= ni) arret("courbe_log: loga = %g; loga2 = %g",loga,loga2) ;
    double loga0 = loga1 + i*dloga ;
    double dv = (Curve_GetYValue(cb)[i+1] - Curve_GetYValue(cb)[i])/dloga ;
    double v  = Curve_GetYValue(cb)[i] + dv*(loga - loga0) ;
    return(v) ;
  }
}



double dcourbe_log(double a,Curve_t* cb)
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double a1 = Curve_GetXRange(cb)[0] ;
  double a2 = Curve_GetXRange(cb)[1] ;
  double loga1 = log10(a1) ;
  double loga2 = log10(a2) ;
  double dloga = (loga2 - loga1)/n_i ;
  double ada = a*pow(10.,dloga) ;
  double da = ada - a ;

  if(a <= a1) return(0.) ; /* pour le cas a = 0 ! */
  else return((courbe_log(a + da,cb) - courbe_log(a - da,cb))*0.5/da) ;
}



double icourbe_log(double a,Curve_t* cb)
/* Return the integral curve computed from cb */ 
{
  int    n_i = Curve_GetNbOfPoints(cb) - 1 ;
  double *x = Curve_GetXRange(cb) ;
  double *y = Curve_GetYValue(cb) ;
  double x1 = x[0] ;
  double x2 = x[1] ;
  double logx1 = log10(x1) ;
  double logx2 = log10(x2) ;
  double dlogx = (logx2 - logx1)/n_i ;
  double loga  = log10(a) ;
  double ln10  = log(10) ;
  double ya = Curve_ComputeValue(cb,a) ;
  double intydx = 0 ;
  int    i ;
  
  for(i = 0 ; logx1 + (i + 1)*dlogx < loga ; i++) {
    double logxi = logx1 + i*dlogx ;
    double xi    = pow(10,logxi) ;
    double logxj = logxi + dlogx ;
    double xj    = pow(10,logxj) ;
    
    intydx += y[i]*xi + y[i + 1]*xj ;
  }
  
  intydx *= dlogx*0.5*ln10 ;
  
  {
    double logxi = logx1 + i*dlogx ;
    double xi    = pow(10,logxi) ;
    
    intydx += (y[i]*xi + ya*a)*(loga - logxi)*0.5*ln10 ;
  }
  
  return(intydx) ;
}

