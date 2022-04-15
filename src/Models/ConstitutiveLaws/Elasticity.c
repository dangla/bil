#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <math.h>

#include "Message.h"
#include "Mry.h"
#include "Math_.h"
#include "Elasticity.h"


static double* (Elasticity_ComputeIsotropicStiffnessTensor)(Elasticity_t*,double*) ;
static double* (Elasticity_ComputeTransverselyIsotropicStiffnessTensor)(Elasticity_t*,double*) ;



Elasticity_t*  (Elasticity_Create)(void)
{
  Elasticity_t* elasty = (Elasticity_t*) Mry_New(Elasticity_t) ;
  
  
  /* Allocation of space for the stiffness tensor */
  {
    double* c = (double*) Mry_New(double[81]) ;
    
    Elasticity_GetStiffnessTensor(elasty) = c ;
  }
  
  /* Allocation of space for the type */
  {
    char* c = (char*) Mry_New(char[Elasticity_MaxLengthOfKeyWord]) ;
    
    Elasticity_GetType(elasty) = c ;
    /* Default = isotropy */
    Elasticity_SetToIsotropy(elasty) ;
  }
  
  /* Allocation of space for the parameters */
  {
    double* c = (double*) Mry_New(double[Elasticity_MaxNbOfParameters]) ;
    
    Elasticity_GetParameter(elasty) = c ;
  }
  
  /* Allocation of space for the stress tensor */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Elasticity_GetStressTensor(elasty) = c ;
  }
  
  return(elasty) ;
}



void  (Elasticity_Delete)(void* self)
{
  Elasticity_t* elasty = (Elasticity_t*) self ;
  
  free(Elasticity_GetStiffnessTensor(elasty)) ;
  free(Elasticity_GetType(elasty)) ;
  free(Elasticity_GetParameter(elasty)) ;
  free(Elasticity_GetStressTensor(elasty)) ;
}



void (Elasticity_SetParameter)(Elasticity_t* elasty,const char* str,double v)
{
  
  if(Elasticity_IsIsotropic(elasty)) {
    if(0) {
    } else if(!strcmp(str,"Young's modulus")) {
      Elasticity_GetYoungModulus(elasty)  = v ;
    } else if(!strcmp(str,"Poisson's ratio")) {
      Elasticity_GetPoissonRatio(elasty)  = v ;
    } else if(!strcmp(str,"bulk modulus")) {
      Elasticity_GetBulkModulus(elasty)   = v ;
    } else if(!strcmp(str,"shear modulus")) {
      Elasticity_GetShearModulus(elasty)  = v ;
    }
    
  } else if(Elasticity_IsTransverselyIsotropic(elasty)) {
    if(0) {
    } else if(!strcmp(str,"Young's modulus")) {
      Elasticity_GetYoungModulus(elasty)  = v ;
    } else if(!strcmp(str,"Poisson's ratio")) {
      Elasticity_GetPoissonRatio(elasty)  = v ;
    } else if(!strcmp(str,"Young's modulus 3")) {
      Elasticity_GetYoungModulus3(elasty) = v ;
    } else if(!strcmp(str,"Poisson's ratio 3")) {
      Elasticity_GetPoissonRatio3(elasty) = v ;
    } else if(!strcmp(str,"shear modulus 3")) {
      Elasticity_GetShearModulus3(elasty) = v ;
    } else if(!strcmp(str,"axis 3")) {
      Elasticity_GetAxis3(elasty)         = v ;
    }
    
  } else {
    Message_RuntimeError("Not known") ;
  }

}



void (Elasticity_SetParameters)(Elasticity_t* elasty,...)
{
  va_list args ;
  
  va_start(args,elasty) ;
  
  if(Elasticity_IsIsotropic(elasty)) {
    double Young   = va_arg(args,double) ;
    double Poisson = va_arg(args,double) ;
    double shear   = Young/(2 + 2*Poisson) ;
    double bulk    = Young/(3 - 6*Poisson) ;
    
    Elasticity_GetYoungModulus(elasty)  = Young ;
    Elasticity_GetPoissonRatio(elasty)  = Poisson ;
    Elasticity_GetBulkModulus(elasty)   = bulk ;
    Elasticity_GetShearModulus(elasty)  = shear ;
    
  } else if(Elasticity_IsTransverselyIsotropic(elasty)) {
    Elasticity_GetYoungModulus(elasty)  = va_arg(args,double) ;
    Elasticity_GetPoissonRatio(elasty)  = va_arg(args,double) ;
    Elasticity_GetYoungModulus3(elasty) = va_arg(args,double) ;
    Elasticity_GetPoissonRatio3(elasty) = va_arg(args,double) ;
    Elasticity_GetShearModulus3(elasty) = va_arg(args,double) ;
    Elasticity_GetAxis3(elasty)         = va_arg(args,double) ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  
  va_end(args) ;
}




double* (Elasticity_ComputeStiffnessTensor)(Elasticity_t* elasty,double* c)
{
  //double* c = Elasticity_GetStiffnessTensor(elasty) ;
  
  if(Elasticity_IsIsotropic(elasty)) {
    return(Elasticity_ComputeIsotropicStiffnessTensor(elasty,c)) ;
  } else if(Elasticity_IsTransverselyIsotropic(elasty)) {
    return(Elasticity_ComputeTransverselyIsotropicStiffnessTensor(elasty,c)) ;
  } else {
    Message_RuntimeError("Not known") ;
  }
  
  return(NULL) ;
}





void (Elasticity_CopyStiffnessTensor)(Elasticity_t* elasty,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* cel = Elasticity_GetStiffnessTensor(elasty) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = cel[i] ;
    }
  }
}




void (Elasticity_PrintStiffnessTensor)(Elasticity_t* elasty)
/** Print the 4th rank elastic tensor.
 **/
{
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  
  printf("\n") ;
  printf("4th rank elastic tensor:\n") ;
  
  Math_PrintStiffnessTensor(c) ;
}




/* Local functions */
double* (Elasticity_ComputeIsotropicStiffnessTensor)(Elasticity_t* elasty,double* c)
/** Compute the 4th rank isotropic elastic tensor in c.
 *  Return c  */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  //double* c = Elasticity_GetStiffnessTensor(elasty) ;
  //double Young   = Elasticity_GetYoungModulus(elasty) ;
  //double Poisson = Elasticity_GetPoissonRatio(elasty) ;
  double bulk    = Elasticity_GetBulkModulus(elasty) ;
  double shear   = Elasticity_GetShearModulus(elasty) ;
  //double twomu   = Young/(1 + Poisson) ;
  //double mu      = (shear > 0) ? shear : Young/(2 + 2*Poisson) ;
  //double k       = (bulk  > 0) ? bulk  : Young/(3 - 6*Poisson) ;
  //double lame    = (Young > 0) ? twomu*Poisson/(1 - 2*Poisson) : k - 2./3.*mu ;
  double mu      = shear ;
  double k       = bulk  ;
  double lame    = k - 2./3.*mu ;
   
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



double* (Elasticity_ComputeTransverselyIsotropicStiffnessTensor)(Elasticity_t* elasty,double* c)
/** Compute the 4th rank transversely isotropic elastic tensor in c.
 *  Inputs are:
 *  axis_3 = direction of orthotropy: 0,1 or 2
 *  Return c  */
{
#define AXIS(I)      (axis_##I)
#define C(i,j,k,l)   (c[((AXIS(i)*3+AXIS(j))*3+AXIS(k))*3+AXIS(l)])
  //double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double Young      = Elasticity_GetYoungModulus(elasty) ;
  double Poisson    = Elasticity_GetPoissonRatio(elasty) ;
  double Young3     = Elasticity_GetYoungModulus3(elasty) ;
  double Poisson3   = Elasticity_GetPoissonRatio3(elasty) ;
  double Shear3     = Elasticity_GetShearModulus3(elasty) ;
  double Axis3      = Elasticity_GetAxis3(elasty) ;
  short int axis_3  = floor(Axis3 + 0.5) ;
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



double* (Elasticity_StiffnessMatrixInVoigtNotation)(Elasticity_t* elasty,double* cv)
/** Compute the stiffness matrix in Voigt's notation from the
 *  4th-order stiffness tensor. This stiffness matrix relates 
 *  the stresses in Voigt's notation: Sij = (S11,S22,S33,S12,S23,S31) to
 *  the strains in Voigt's notation:  Eij = (E11,E22,E33,2E12,2E23,2E31)
 *  On inputs:
 *  - elasty should contain the initial 4th-order stiffness tensor,
 *  - cv should point to an allocated space of 36 doubles.
 *  On output:
 *  - cv will contain the stiffness matrix in Voigt's notation
 **/
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
#define CV(i,j)     (cv[(i)*6+(j)])
/*
 *                  |0 3 5|
 * VI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define VI(i,j)     (voigtindex[(i)*3 + (j)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  int voigtindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ij = VI(i,j) ;
      int k ;
      
      for(k = 0 ; k < 3 ; k++) {
        int l  = (k + 1) % 3 ;
        int kl = VI(k,l) ;
        
        CV(i,k)   = C(i,i,k,k) ;
        CV(i,kl)  = 0.5  * (C(i,i,k,l) + C(i,i,l,k)) ;
        CV(kl,i)  = 0.5  * (C(k,l,i,i) + C(l,k,i,i)) ;
        CV(ij,kl) = 0.25 * (C(i,j,k,l) + C(j,i,k,l) + C(i,j,l,k) + C(j,i,l,k)) ;
      }
    }
  }
  
  return(cv) ;
  
#undef VI
#undef CV
#undef C
}





double (Elasticity_ComputeElasticEnergy)(Elasticity_t* elasty,const double* strain)
/** Return the elastic energy associated to the strain tensor.
 **/
{
  double energy = 0 ;
  
  #define C(i,j)  (c[(i)*9+(j)])
  {
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    int i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double stress = 0 ;
      int j ;
      
      for(j = 0 ; j < 9 ; j++) {
        stress += C(i,j) * strain[j] ;
      }
      
      energy += strain[i] * stress ;
    }
    
    energy *= 0.5 ;
  }
  #undef C
  
  return(energy) ;
}





double* (Elasticity_ComputeStressTensor)(Elasticity_t* elasty,const double* strain,double* stress0)
/** Return the elastic stress tensor associated to the strain tensor.
 **/
{
  #define C(i,j)  (c[(i)*9+(j)])
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  double* stress = (stress0) ? stress0 : Elasticity_GetStressTensor(elasty) ;
  int i ;
    
  for(i = 0 ; i < 9 ; i++) {
    int j ;
      
    stress[i] = 0 ;
      
    for(j = 0 ; j < 9 ; j++) {
      stress[i] += C(i,j) * strain[j] ;
    }
  }
  #undef C
  
  return(stress) ;
}

