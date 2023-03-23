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
  
  /* Allocation of space for the compliance tensor */
  {
    double* c = (double*) Mry_New(double[81]) ;
    
    Elasticity_GetComplianceTensor(elasty) = c ;
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
  
  {
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Elasticity_GetComplianceTensor(elasty) ;
    
    if(c) free(c) ;
  }
  
  {
    char* type = Elasticity_GetType(elasty) ;
    
    if(type) free(type) ;
  }
  
  {
    double* c = Elasticity_GetParameter(elasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Elasticity_GetStressTensor(elasty) ;
    
    if(c) free(c) ;
  }
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
    
    if(Young < 0) {
      Message_RuntimeError("Elasticity_SetParameters: negative Young's modulus") ;
    }
    
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
    Message_RuntimeError("Elasticity_SetParameters: illegal elasticity") ;
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





double* (Elasticity_CopyStiffnessTensor)(Elasticity_t* elasty,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* cel = Elasticity_GetStiffnessTensor(elasty) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = cel[i] ;
    }
  }
  
  return(c) ;
}





double* (Elasticity_CopyComplianceTensor)(Elasticity_t* elasty,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* cel = Elasticity_GetComplianceTensor(elasty) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = cel[i] ;
    }
  }
  
  return(c) ;
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



double* (Elasticity_ConvertStiffnessMatrixInto6x6Matrix)(double* c)
/** Convert the 4th-order stiffness tensor into a 6x6 matrix using
 *  the Mandel's notation. This 6x6 Mandel's matrix relates 
 *  the stresses in Mandel's notation: Sij = (S11,S22,S33,sr2*S12,sr2*S23,sr2*S31) to
 *  the strains in Mandel's notation:  Eij = (E11,E22,E33,sr2*E12,sr2*E23,sr2*E31)
 *  where sr2 stands for the square root of 2.
 *  On inputs:
 *  - c should point to an allocated space of 81 doubles containing
 *    the 81 components of the 4th-order stiffness matrix.
 *  On output:
 *  - c is replaced by the 2nd-order stiffness matrix in Mandel's notation.
 *    Only the first 36 doubles are modified. 
 *    More precisely the stiffness matrix in Mandel's notation writes
 *         C1111      C1122      C1133  sr2*C1112  sr2*C1123  sr2*C1131
 *         C2211      C2222      C2233  sr2*C2212  sr2*C2223  sr2*C2231
 *         C3311      C3322      C3333  sr2*C3312  sr2*C3323  sr2*C3331
 *     sr2*C1211  sr2*C1222  sr2*C1233    2*C1212    2*C1223    2*C1231
 *     sr2*C2311  sr2*C2322  sr2*C2333    2*C2312    2*C2323    2*C2331
 *     sr2*C3111  sr2*C3122  sr2*C3133    2*C3112    2*C3123    2*C3131
 *  ref: 
 *  T. Manik, A natural vector/matrix notation applied in an efficient 
 *  and robust return-mapping algorithm for advanced yield functions, 
 *  European Journal of Mechanics / A Solids, 90 (2021).
 **/
{
#define C4(i,j,k,l) (c[(((i)*3+(j))*3+(k))*3+(l)])
#define C2(i,j)     (c2[(i)*6+(j)])
/*
 *                  |0 3 5|
 * MI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define MI(i,j)     (mapindex[(i)*3 + (j)])
  double  c2[36] ;
  double z5sr2 = 0.5*sqrt(2) ; /* = 1/sqrt(2) */
  double z5    = 0.5 ;
  int mapindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ii = MI(i,i) ;
      int ij = MI(i,j) ;
      int k ;
      
      for(k = 0 ; k < 3 ; k++) {
        int l  = (k + 1) % 3 ;
        int kk = MI(k,k) ;
        int kl = MI(k,l) ;
        
        C2(ii,kk) = C4(i,i,k,k) ;
        C2(ii,kl) = z5sr2 * (C4(i,i,k,l) + C4(i,i,l,k)) ;
        C2(kl,ii) = z5sr2 * (C4(k,l,i,i) + C4(l,k,i,i)) ;
        C2(ij,kl) = z5 * (C4(i,j,k,l) + C4(j,i,k,l) + C4(i,j,l,k) + C4(j,i,l,k)) ;
      }
    }
  }
  
  {
    int    i ;

    for(i = 0 ; i < 36 ; i++) {
      c[i] = c2[i] ;
    }
  }
  
  return(c) ;
  
#undef MI
#undef C2
#undef C4
}



double* (Elasticity_Convert6x6MatrixIntoStiffnessMatrix)(double* c)
/** Convert a 6x6 stiffness matrix in Mandel's notation into a 
 *  4th-order stiffness matrix. This 6x6 Mandel's matrix relates 
 *  the stresses in Mandel's notation: Sij = (S11,S22,S33,sr2*S12,sr2*S23,sr2*S31) to
 *  the strains in Mandel's notation:  Eij = (E11,E22,E33,sr2*E12,sr2*E23,sr2*E31)
 *  where sr2 stands for the square root of 2.
 *  On inputs:
 *  - c should point to an allocated space of 81 doubles. 
 *    The first 36 doubles should contain the matrix components 
 *    in Mandel's notation
 *  On output:
 *  - c is replaced by the 4th-order symmetric matrix with 81 components.
 *    More precisely the stiffness matrix in Mandel's notation writes
 *         C1111      C1122      C1133  sr2*C1112  sr2*C1123  sr2*C1131
 *         C2211      C2222      C2233  sr2*C2212  sr2*C2223  sr2*C2231
 *         C3311      C3322      C3333  sr2*C3312  sr2*C3323  sr2*C3331
 *     sr2*C1211  sr2*C1222  sr2*C1233    2*C1212    2*C1223    2*C1231
 *     sr2*C2311  sr2*C2322  sr2*C2333    2*C2312    2*C2323    2*C2331
 *     sr2*C3111  sr2*C3122  sr2*C3133    2*C3112    2*C3123    2*C3131
 *  ref: 
 *  T. Manik, A natural vector/matrix notation applied in an efficient 
 *  and robust return-mapping algorithm for advanced yield functions, 
 *  European Journal of Mechanics / A Solids, 90 (2021).
 **/
{
#define C4(i,j,k,l) (c4[(((i)*3+(j))*3+(k))*3+(l)])
#define C2(i,j)     (c[(i)*6+(j)])
/*
 *                  |0 3 5|
 * MI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define MI(i,j)     (mapindex[(i)*3 + (j)])
  double c4[81] ;
  double z5sr2 = 0.5*sqrt(2) ; /* = 1/sqrt(2) */
  double z5    = 0.5 ;
  int mapindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ii = MI(i,i) ;
      int ij = MI(i,j) ;
      int k ;
      
      for(k = 0 ; k < 3 ; k++) {
        int l  = (k + 1) % 3 ;
        int kk = MI(k,k) ;
        int kl = MI(k,l) ;
        
        C4(i,i,k,k) = C2(ii,kk) ;
        C4(i,i,k,l) = z5sr2 * C2(ii,kl) ;
        C4(i,i,l,k) = z5sr2 * C2(ii,kl) ;
        C4(k,l,i,i) = z5sr2 * C2(kl,ii) ;
        C4(l,k,i,i) = z5sr2 * C2(kl,ii) ;
        C4(i,j,k,l) = z5 * C2(ij,kl) ;
        C4(j,i,k,l) = z5 * C2(ij,kl) ;
        C4(i,j,l,k) = z5 * C2(ij,kl) ;
        C4(j,i,l,k) = z5 * C2(ij,kl) ;
      }
    }
  }
  
  {
    int    i ;

    for(i = 0 ; i < 81 ; i++) {
      c[i] = c4[i] ;
    }
  }
  
  return(c) ;
  
#undef MI
#undef C2
#undef C4
}



double* (Elasticity_ConvertStressTensorInto6TermStressVector)(double* c)
/** Convert the 2nd-order symmetric stress tensor into a 6x1 column 
 *  matrix using the Mandel's notation. This 6x1 Mandel's column matrix 
 *  relates the stress tensor to the stresses in Mandel's notation: 
 *  Sij -> (S11,S22,S33,sr2*S12,sr2*S23,sr2*S31)
 *  where sr2 stands for the square root of 2.
 *  On inputs:
 *  - c should point to an allocated space of 9 doubles containing
 *    the 9 components of the 2nd-order symmetric stress tensor.
 *  On output:
 *  - c is replaced by the column stress vector in Mandel's notation.
 *    Only the first 6 doubles are modified.
 **/
{
#define C2(i,j)     (c[(i)*3+(j)])
#define C1(i)       (c1[i])
/*
 *                  |0 3 5|
 * MI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define MI(i,j)     (mapindex[(i)*3 + (j)])
  double  c1[6] ;
  double z5sr2 = 0.5*sqrt(2) ;
  int mapindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ii = MI(i,i) ;
      int ij = MI(i,j) ;
      
      C1(ii) = C2(i,i) ;
      C1(ij) = z5sr2 * (C2(i,j) + C2(j,i)) ;
    }
  }
  
  {
    int    i ;

    for(i = 0 ; i < 6 ; i++) {
      c[i] = c1[i] ;
    }
  }
  
  return(c) ;
  
#undef MI
#undef C2
#undef C1
}



double* (Elasticity_Convert6TermStressVectorIntoStressTensor)(double* c)
/** Convert a 6 term column matrix representing a symmetric vector 
 *  using the Mandel's notation into a 2nd-order symmetric stress tensor. 
 *  The Mandel's notation relates the stress tensor to the 6 term vector
 *  of stresses: 
 *  Sij -> (S11,S22,S33,sr2*S12,sr2*S23,sr2*S31)
 *  where sr2 stands for the square root of 2.
 *  On input:
 *  - c should point to an allocated space of 9 doubles containing in
 *    the first 6 terms the column vector in Mandel's notation.
 *  On output:
 *  - c is replaced by the 9 components of the 2nd-order symmetric
 *    stress tensor.
 **/
{
#define C2(i,j)     (c2[(i)*3+(j)])
#define C1(i)       (c[i])
/*
 *                  |0 3 5|
 * MI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define MI(i,j)     (mapindex[(i)*3 + (j)])
  double  c2[9] ;
  double z5sr2 = 0.5*sqrt(2) ;
  int mapindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ii = MI(i,i) ;
      int ij = MI(i,j) ;
      
      C2(i,i) = C1(ii) ;
      C2(i,j) = z5sr2 * C1(ij) ;
      C2(j,i) = z5sr2 * C1(ij) ;
    }
  }
  
  {
    int    i ;

    for(i = 0 ; i < 9 ; i++) {
      c[i] = c2[i] ;
    }
  }
  
  return(c) ;
  
#undef MI
#undef C2
#undef C1
}



double* (Elasticity_InvertStiffnessMatrix)(double* c)
/** Replace a 4th-order stiffness tensor by its inverse.
 **/
{
  
  Elasticity_ConvertStiffnessMatrixInto6x6Matrix(c) ;

  Math_InvertMatrix(c,6) ;
  
  Elasticity_Convert6x6MatrixIntoStiffnessMatrix(c) ;
    
  return(c) ;
}



#if 0
double* (Elasticit_IsotropicStiffnessTensor)(const double k, const double g,double* c)
/** Compute the 4th rank isotropic elastic tensor from k and g.
 *  Return c  */
{
#define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
  double lame    = k - 2./3.*g ;
   
  {
    int    i ;

    for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
    for(i = 0 ; i < 3 ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) {
        C(i,i,j,j) += lame ;
        C(i,j,i,j) += g ;
        C(i,j,j,i) += g ;
      }
    }
  }
  
  return(c) ;
#undef C
}
#endif



#if 0
double* (Elasticity_StiffnessMatrixInVoigtNotation)(double* c)
/** Compute the stiffness matrix in Voigt's notation from the
 *  4th-order stiffness tensor. This stiffness Voigt's matrix relates 
 *  the stresses in Voigt's notation: Sij = (S11,S22,S33,S12,S23,S31) to
 *  the strains in Voigt's notation:  Eij = (E11,E22,E33,2E12,2E23,2E31)
 *  On inputs:
 *  - c should point to an allocated space of 81 doubles containing
 *    the 81 components of the 4th-order stiffness matrix.
 *  On output:
 *  - c is replaced by the 2nd-order stiffness matrix in Voigt's notation.
 *    Only the first 36 doubles are modified.
 **/
{
#define C4(i,j,k,l) (c[(((i)*3+(j))*3+(k))*3+(l)])
#define C2(i,j)     (c2[(i)*6+(j)])
/*
 *                  |0 3 5|
 * MI(i,j) maps to  |3 1 4|
 *                  |5 4 2|
 */
#define MI(i,j)     (mapindex[(i)*3 + (j)])
  double  c2[36] ;
  int mapindex[9] = {0,3,5,3,1,4,5,4,2} ;
  
  {
    int    i ;

    for(i = 0 ; i < 3 ; i++) {
      int j  = (i + 1) % 3 ;
      int ii = MI(i,i) ;
      int ij = MI(i,j) ;
      int k ;
      
      for(k = 0 ; k < 3 ; k++) {
        int l  = (k + 1) % 3 ;
        int kk = MI(k,k) ;
        int kl = MI(k,l) ;
        
        C2(ii,kk) = C4(i,i,k,k) ;
        C2(ii,kl) = 0.5  * (C4(i,i,k,l) + C4(i,i,l,k)) ;
        C2(kl,ii) = 0.5  * (C4(k,l,i,i) + C4(l,k,i,i)) ;
        C2(ij,kl) = 0.25 * (C4(i,j,k,l) + C4(j,i,k,l) + C4(i,j,l,k) + C4(j,i,l,k)) ;
      }
    }
  }
  
  {
    int    i ;

    for(i = 0 ; i < 36 ; i++) {
      c[i] = c2[i] ;
    }
  }
  
  return(c) ;
  
#undef MI
#undef C2
#undef C4
}
#endif



#if 0
static int Elasticity_Test(int,char**) ;

#include <time.h>
#include "Math_.h"

int Elasticity_Test(int argc,char** argv)
{
  #define  A(i,j)  (a[(i)*3+(j)])
  double bulk  = (argc > 1) ? (double) atof(argv[1]) : 1 ;
  double shear = (argc > 2) ? (double) atof(argv[2]) : 1 ;
  double poisson = (3*bulk - 2*shear)/(6*bulk + 2*shear) ;
  double young = 9*bulk*shear/(3*bulk+shear) ;
  Elasticity_t* elasty = Elasticity_Create() ;
  
  #if 0
  printf("Bulk's modulus: %g\n",bulk) ;
  printf("Shear modulus: %g\n",shear) ;
  printf("Young's modulus: %g\n",young) ;
  printf("Poisson's ratio: %g\n",poisson) ;
  
  Elasticity_SetToIsotropy(elasty) ;
  Elasticity_SetParameters(elasty,young,poisson) ;
    
  Elasticity_ComputeStiffnessTensor(elasty,c) ;
  #endif
  
  {
    double c[81] ;
    double invc[81] ;
    double id[81] ;
    
    /* Fill a 6x6 matrix randomly */
    {
      int i ;
      int rmax = RAND_MAX / 2 ;
      
      srand(time(NULL)) ;
      srand(rand()) ;
      
      #define C2(i,j)  c[(i)*6 + (j)]
      for(i = 0 ; i < 6 ; i++) {
        int j ;
        
        C2(i,i) = ((double) (rand() - rmax))/rmax ;
        
        for(j = 0 ; j < 6 ; j++) {
          C2(i,j) = ((double) (rand() - rmax))/rmax ;
        }
      }
      #undef C2
    }
    
    /* Convert it into a 4th-order stiffness matrix */
    Elasticity_Convert6x6MatrixIntoStiffnessMatrix(c) ;
    
    /* Save it in invc */
    {
      int i ;
      
      for(i = 0 ; i < 81 ; i++) {
        invc[i] = c[i] ;
      }
    }
  
    printf("original c:\n") ;
    Math_PrintStiffnessTensor(c) ;
    
    Elasticity_InvertStiffnessMatrix(invc) ;
    
    printf("inverse c:\n") ;
    Math_PrintStiffnessTensor(invc) ;
    
    /* check */
    {
      int i ;
      
      #define ID(i,j)    id[(i)*9 + (j)]
      #define C(i,j)     c[(i)*9 + (j)]
      #define INVC(i,j)  invc[(i)*9 + (j)]
      for(i = 0 ; i < 9 ; i++) {
        int j ;
        
        for(j = 0 ; j < 9 ; j++) {
          int k ;
          
          ID(i,j) = 0  ;
          
          for(k = 0 ; k < 9 ; k++) {
            ID(i,j) += C(i,k)*INVC(k,j)  ;
          }
        }
      }
      #undef ID
      #undef C
      #undef INVC
    
      printf("check original * inverse?:\n") ;
      Math_PrintStiffnessTensor(id) ;
    }
    
    #if 0
    {
      bulk  = 1/(9*bulk) ;
      shear = 1/(4*shear) ;
      
      poisson = (3*bulk - 2*shear)/(6*bulk + 2*shear) ;
      young = 9*bulk*shear/(3*bulk+shear) ;
      
      Elasticity_SetParameters(elasty,young,poisson) ;
      Elasticity_ComputeIsotropicStiffnessTensor(elasty,c) ;
      
      printf("check c\n") ;
      
      Math_PrintStiffnessTensor(c) ;
    }
    #endif
  }
  
  return(0) ;
}


/*
 * Compilation: 
 * g++ -gdwarf-2 -g3  -L/home/dangla/Documents/Softwares/bil/bil-master/lib -Wl,-rpath=/home/dangla/Documents/Softwares/bil/bil-master/lib -lbil-2.8.8-Debug -o out -lgfortran
*/
int main(int argc, char** argv)
{
  Session_Open() ;
  
  Elasticity_Test(argc,argv) ;
  
  Session_Close() ;
  return(0) ;
}
#endif
