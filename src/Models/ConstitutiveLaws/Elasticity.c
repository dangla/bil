#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "Elasticity.h"


static double* Elasticity_ComputeIsotropicStiffnessTensor(Elasticity_t* elasty,const double,const double) ;
static double* Elasticity_ComputeTransverselyIsotropicStiffnessTensor(Elasticity_t* elasty,const double,const double,const double,const double,const double,const short int) ;



Elasticity_t*  (Elasticity_Create)(void)
{
  Elasticity_t* elasty = (Elasticity_t*) malloc(sizeof(Elasticity_t)) ;
  
  assert(elasty) ;
  
  /* Allocation of space for the stiffness tensor */
  {
    size_t sz = 81*sizeof(double) ;
    double* c = (double*) malloc(sz) ;
    
    assert(c) ;
    
    Elasticity_GetStiffnessTensor(elasty) = c ;
  }
  
  return(elasty) ;
}



void  (Elasticity_Delete)(void* self)
{
  Elasticity_t** pelasty = (Elasticity_t**) self ;
  Elasticity_t*  elasty  = *pelasty ;
  
  {
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    free(c) ;
  }
  
  free(*pelasty) ;
  *pelasty = NULL ;
}




double* Elasticity_ComputeIsotropicStiffnessTensor(Elasticity_t* elasty,const double Young,const double Poisson)
/** Compute the 4th rank isotropic elastic tensor in c.
 *  Return c  */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double twomu   = Young/(1 + Poisson) ;
  double mu      = 0.5*twomu ;
  double lame    = twomu*Poisson/(1 - 2*Poisson) ;
   
  {
    int    i ;

    for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
    for(i = 0 ; i < 3 ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) {
        C(i,i,j,j) += lame ;
        C(i,j,i,j) += mu ;
        C(i,j,j,i) += mu ;
      }
    }
  }
  
  return(c) ;
#undef C
}



double Elasticity_ComputeYoungModulus(Elasticity_t* elasty)
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double mu     = C(0,1,0,1) ;
  double lame   = C(0,0,0,0) - 2 * mu ;
  double poisson = 0.5 * lame / (mu + lame) ;
  double young  = 2 * mu * (1 + poisson) ;
  
  return(young) ;
#undef C
}



double Elasticity_ComputePoissonRatio(Elasticity_t* elasty)
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double mu     = C(0,1,0,1) ;
  double lame   = C(0,0,0,0) - 2 * mu ;
  double poisson = 0.5 * lame / (mu + lame) ;
  
  return(poisson) ;
#undef C
}



double Elasticity_ComputeShearModulus(Elasticity_t* elasty)
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double mu = C(0,1,0,1) ;
  
  return(mu) ;
#undef C
}



double* Elasticity_ComputeTransverselyIsotropicStiffnessTensor(Elasticity_t* elasty,const double Young,const double Young3,const double Poisson,const double Poisson3,const double Shear3,const short int axis_3)
/** Compute the 4th rank transversely isotropic elastic tensor in c.
 *  Inputs are:
 *  axis_3 = direction of orthotropy: 0,1 or 2
 *  Return c  */
{
#define AXIS(I)      (axis_##I)
#define C(i,j,k,l)   (c[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  short int axis_1  = (axis_3 + 1) % 3 ;
  short int axis_2  = (axis_3 + 2) % 3 ;
  double twomu1     = Young/(1 + Poisson) ;
  double twomu3     = 2*Shear3 ;
  double mu1        = 0.5*twomu1 ;
  double mu3        = 0.5*twomu3 ;
  double lamu1      = (1 - Poisson)/Young - 2*Poisson3*Poisson3/Young3 ;
  double lam1       = 0.5/lamu1 - 0.5*twomu1 ;
  double lam2       = Poisson3/lamu1 ;
  double lam3       = Young3 + 2*Poisson3*lam2 - twomu3 ;
   
  {
    int    i ;

    for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
  }
  
  {
    C(1,1,1,1) = lam1 + twomu1 ;
    C(1,1,2,2) = lam1 ;
    C(1,1,3,3) = lam2 ;
    
    C(2,2,1,1) = lam1 ;
    C(2,2,2,2) = lam1 + twomu1 ;
    C(2,2,3,3) = lam2 ;
    
    C(3,3,1,1) = lam2 ;
    C(3,3,2,2) = lam2 ;
    C(3,3,3,3) = lam3 + twomu3 ;
    
    C(1,2,1,2) = mu1 ;
    C(1,2,2,1) = mu1 ;
    C(2,1,1,2) = mu1 ;
    C(2,1,2,1) = mu1 ;
    
    C(1,3,1,3) = mu3 ;
    C(1,3,3,1) = mu3 ;
    C(3,1,1,3) = mu3 ;
    C(3,1,3,1) = mu3 ;
    
    C(2,3,2,3) = mu3 ;
    C(2,3,3,2) = mu3 ;
    C(3,2,2,3) = mu3 ;
    C(3,2,3,2) = mu3 ;
  }
  
  return(c) ;
  
#undef C
#undef AXIS
}



#if 0
void Elasticity_PrintStiffnessTensor(Elasticity_t* elasty)
/** Print the 4th rank elastic tensor.
 **/
{
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  
  printf("\n") ;
  printf("4th rank elastic tensor:\n") ;
  
  {
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      int j = i - (i/3)*3 ;
        
      printf("C%d%d--:",i/3 + 1,j + 1) ;
        
      for (j = 0 ; j < 9 ; j++) {
        printf(" % e",c[i*9 + j]) ;
      }
        
      printf("\n") ;
    }
  }
}
#endif
