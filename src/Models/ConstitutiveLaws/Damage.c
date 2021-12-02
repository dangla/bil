#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Mry.h"
#include "Message.h"
#include "DataFile.h"
#include "Tools/Math.h"
#include "Damage.h"


/* Mazars */
static Damage_ComputeFunctionGradients_t    Damage_ComputeFunctionGradientsMazars ;
static Damage_ReturnMapping_t               Damage_ReturnMappingMazars ;
/* Marigo-Jirasek */
static Damage_ComputeFunctionGradients_t    Damage_ComputeFunctionGradientsMarigoJirasek ;
static Damage_ReturnMapping_t               Damage_ReturnMappingMarigoJirasek ;





Damage_t*  (Damage_Create)(void)
{
  Damage_t* damage = (Damage_t*) Mry_New(Damage_t) ;
  
  
  /* Allocation of space for the code name of the model */
  {
    char* name = (char*) Mry_New(char[Damage_MaxLengthOfKeyWord]) ;
    
    Damage_GetCodeNameOfModel(damage) = name ;
  }
  
  /* Allocation of space for the yield function gradient */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Damage_GetYieldFunctionGradient(damage) = c ;
  }
  
  /* Allocation of space for the potential function gradient */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Damage_GetPotentialFunctionGradient(damage) = c ;
  }
  
  /* Allocation of space for the hardening variable */
  {
    double* c = (double*) Mry_New(double[Damage_MaxNbOfHardeningVariables]) ;
    
    Damage_GetHardeningVariable(damage) = c ;
  }
  
  /* Allocation of space for the hardening modulus */
  {
    double* c = (double*) Mry_New(double) ;
    
    Damage_GetHardeningModulus(damage) = c ;
  }
  
  /* Allocation of space for Fji*Cijkl */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Damage_GetFjiCijkl(damage) = c ;
  }
  
  /* Allocation of space for Cijkl*Glk */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Damage_GetCijklGlk(damage) = c ;
  }
  
  /* Allocation of space for the tangent stiffness tensor */
  {
    double* c = (double*) Mry_New(double[81]) ;
    
    Damage_GetStiffnessTensor(damage) = c ;
  }
  
  /* Elasticity */
  Damage_GetElasticity(damage) = Elasticity_Create() ;
  
  /* Allocation of space for the parameters */
  {
    double* c = (double*) Mry_New(double[Damage_MaxNbOfParameters]) ;
    
    Damage_GetParameter(damage) = c ;
  }
  
  return(damage) ;
}



void  (Damage_Delete)(void* self)
{
  Damage_t* damage = (Damage_t*) self ;
  
  {
    char* name = Damage_GetCodeNameOfModel(damage) ;
    free(name) ;
  }
  
  {
    double* c = Damage_GetYieldFunctionGradient(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetPotentialFunctionGradient(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetHardeningVariable(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetHardeningModulus(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetFjiCijkl(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetCijklGlk(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetStiffnessTensor(damage) ;
    free(c) ;
  }
  
  {
    Elasticity_t* elasty = Damage_GetElasticity(damage) ;
    
    Elasticity_Delete(elasty) ;
    free(elasty) ;
  }
  
  {
    double* c = Damage_GetParameter(damage) ;
    free(c) ;
  }
}



void Damage_Initialize(Damage_t* damage)
{
  
  if(Damage_IsMazars(damage)) {
    Damage_GetComputeFunctionGradients(damage) = Damage_ComputeFunctionGradientsMazars ;
    Damage_GetReturnMapping(damage)            = Damage_ReturnMappingMazars ;
      
  } else if(Damage_IsMarigoJirasek(damage)) {
    Damage_GetComputeFunctionGradients(damage) = Damage_ComputeFunctionGradientsMarigoJirasek ;
    Damage_GetReturnMapping(damage)            = Damage_ReturnMappingMarigoJirasek ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  
}



void Damage_SetParameter(Damage_t* damage,const char* str,double v)
{
  
  if(Damage_IsMazars(damage)) {
    if(0) {
    } else if(!strcmp(str,"max elastic strain")) {
      Damage_GetStrainAtUniaxialTensileStrength(damage)  = v ;
      Damage_GetHardeningVariable(damage)[0]  = v ;
    } else if(!strcmp(str,"A_c")) {
      Damage_GetA_c(damage)  = v ;
    } else if(!strcmp(str,"A_t")) {
      Damage_GetA_t(damage)  = v ;
    } else if(!strcmp(str,"B_c")) {
      Damage_GetB_c(damage)  = v ;
    } else if(!strcmp(str,"B_t")) {
      Damage_GetB_t(damage)  = v ;
    }
    
  } else if(Damage_IsMarigoJirasek(damage)) {
    Elasticity_t* elasty = Damage_GetElasticity(damage) ;
    double young = Elasticity_GetYoungModulus(elasty) ;
    
    if(0) {
    } else if(!strcmp(str,"uniaxial tensile strength")) {
      double ft = v ;
      double strain0 = ft/young ;
      double g0   = 0.5*young*strain0*strain0 ;
      
      Damage_GetUniaxialTensileStrength(damage)  = ft ;
      Damage_GetHardeningVariable(damage)[0] = g0 ;
    } else if(!strcmp(str,"fracture energy")) {
      Damage_GetFractureEnergy(damage)  = v ;
    } else if(!strcmp(str,"crack band width")) {
      Damage_GetCrackBandWidth(damage)  = v ;
    }
    
  } else {
    Message_RuntimeError("Not known") ;
  }

}



void Damage_SetParameters(Damage_t* damage,...)
{
  va_list args ;
  
  va_start(args,damage) ;
  
  if(Damage_IsMazars(damage)) {
    double strain0  = va_arg(args,double) ;
    double A_c  = va_arg(args,double) ;
    double A_t  = va_arg(args,double) ;
    double B_c  = va_arg(args,double) ;
    double B_t  = va_arg(args,double) ;
    
    Damage_GetStrainAtUniaxialTensileStrength(damage) = strain0 ;
    Damage_GetA_c(damage)  = A_c ;
    Damage_GetA_t(damage)  = A_t ;
    Damage_GetB_c(damage)  = B_c ;
    Damage_GetB_t(damage)  = B_t ;
    
    Damage_GetHardeningVariable(damage)[0] = strain0 ;
    
  } else if(Damage_IsMarigoJirasek(damage)) {
    double ft  = va_arg(args,double) ;
    double Gf  = va_arg(args,double) ;
    double w   = va_arg(args,double) ;
    Elasticity_t* elasty = Damage_GetElasticity(damage) ;
    double E       = Elasticity_GetYoungModulus(elasty) ;
    double strain0 = ft/E ;
    double g0      = 0.5*E*strain0*strain0 ;
    
    Damage_GetUniaxialTensileStrength(damage)   = ft ;
    Damage_GetFractureEnergy(damage)            = Gf ;
    Damage_GetCrackBandWidth(damage)            = w ;
    
    Damage_GetHardeningVariable(damage)[0]      = g0 ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  

  va_end(args) ;
}





void Damage_CopyStiffnessTensor(Damage_t* damage,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* cel = Damage_GetStiffnessTensor(damage) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = cel[i] ;
    }
  }
}









double Damage_UpdateTangentStiffnessTensor(Damage_t* damage,double* c)
/** Update the 4th rank tangent damage tensor in c.
 *  On input c should point to the elastic stiffness tensor. 
 *  On output c is updated to the tangent damage stiffness tensor.
 *  Other inputs, included in damage, are: 
 *  dfsde is the yield function gradient
 *  dgsde is the potential function gradient
 *  hm    is the hardening modulus
 *  Tensor c is then updated as 
 *  C(i,j,k,l) = C^el(i,j,k,l) - dgsde(i,j) * dfsde(k,l) / hm
 *  Return the inverse of hm: 1/hm */
{
#define CD(i,j)   ((c)[(i)*9+(j)])
#define CEL(i,j)  ((c)[(i)*9+(j)])
//#define CEL(i,j)  ((cel)[(i)*9+(j)])
  //Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  //double* cel    = Elasticity_GetStiffnessTensor(elasty) ;
  double* dfsde  = Damage_GetYieldFunctionGradient(damage) ;
  double* dgsde  = Damage_GetPotentialFunctionGradient(damage) ;
  double* hm     = Damage_GetHardeningModulus(damage) ;
  double  det ;
  
  /* Tangent stiffness tensor */
  {
      int i ;
      
      det = hm[0] ;
        
      if(det > 0.) {
        det = 1./det ;
      } else {
            
        printf("\n") ;
        printf("dF = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dfsde[i]) ;
        }
        printf("\n") ;
        printf("dG = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dgsde[i]) ;
        }
        printf("\n") ;
        printf("hm = %e\n",det) ;
        printf("\n") ;
        
        arret("Damage_UpdateTangentStiffnessTensor") ;
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          CD(i,j) = CEL(i,j) - dgsde[i]*dfsde[j]*det ;
        }
      }
  }
  
  return(det) ;
#undef CD
#undef CEL
}




void Damage_PrintStiffnessTensor(Damage_t* damage)
/** Print the 4th rank tangent damage tensor.
 **/
{
  double* c = Damage_GetStiffnessTensor(damage) ;
  
  printf("\n") ;
  printf("4th rank tangent damage tensor:\n") ;
  
  Math_PrintStiffnessTensor(c) ;
}







/* Local functions */
double Damage_ComputeFunctionGradientsMazars(Damage_t* damage,const double* strain,const double* d,const double* hardv)
/** Mazars criterion.
 *  J. Mazars, A description of micro- and macroscale damage of concrete structures,
 *  Engineering Fracture Mechanics (1986), 25(5/6): 729-737.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest equivalent strain (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the Poisson's ratio (poisson),
 *  the limit elastic strain at the peak of uniaxial tensile stress test,
 *  the parameters of Mazars' model: A_t,B_t,A_c,B_c.
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  dgsde = derivative of the potential function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double* dfsde  = Damage_GetYieldFunctionGradient(damage) ;
  double* dgsde  = Damage_GetPotentialFunctionGradient(damage) ;
  double* hm     = Damage_GetHardeningModulus(damage) ;
  double kappa0  = Damage_GetStrainAtUniaxialTensileStrength(damage) ;
  double A_t     = Damage_GetA_t(damage) ;
  double B_t     = Damage_GetB_t(damage) ;
  double A_c     = Damage_GetA_c(damage) ;
  double B_c     = Damage_GetB_c(damage) ;

  double* strain_val ; //= Math_ComputePrincipalStresses(strain) ;
  double strain_dir[9] ;
  double strain_pos[3] ;
  double strain_eq ;
  double crit ;
  
  /* Principal values and principal directions */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      strain_dir[i] = strain[i] ;
    }
    
    strain_val = Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix(strain_dir,'r') ;
  }
  
  /*
    Positive principal strains
  */
  {
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_pos[i] = Math_Max(strain_val[i],0) ;
    }
  }
  
  /*
    Equivalent tensile strain
  */
  {
    int    i ;
    
    strain_eq = 0 ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_eq += strain_pos[i]*strain_pos[i] ;
    }
    
    strain_eq = sqrt(strain_eq) ;
  }
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    
    crit = strain_eq - kappa ;
  }
  
  /*
    Function gradients
    For an isotropic function of the strain tensor both
    the strain and the function gradient have the same eigenvectors.
    Denote with (e1,e2,e3) the eigenvalues of the strain tensor
    and consider the function f(e1,e2,e3): the eigenvalues
    of the function gradient are df/de1, df/de2, df/de3.
  */
  {
    int    i ;

    #define DF(i,j)   dfsde[3*(i) + (j)]
    for(i = 0 ; i < 3 ; i++) {
      int j ;
      
      for(j = 0 ; j < 3 ; j++) {
        int k ;
        
        DF(i,j) = 0 ;
        
        for(k = 0 ; k < 3 ; k++) {
          double dij = strain_dir[3*(k) + i] * strain_dir[3*(k) + j] ;
          
          DF(i,j) += strain_pos[k] / strain_eq * dij ;
        }
      }
    }
    #undef DF
    
    for(i = 0 ; i < 9 ; i++) {
      dgsde[i] = dfsde[i] ;
    }
  }
  
  /* Hardening modulus: it is the derivative of kappa w.r.t. damage */
  {
    double kappa = hardv[0] ;
    double kappa2 = kappa * kappa ;
    double dd_t = kappa0 / kappa2 * (1 - A_t) + A_t * B_t * exp(-B_t * (kappa - kappa0)) ;
    double dd_c = kappa0 / kappa2 * (1 - A_c) + A_c * B_c * exp(-B_c * (kappa - kappa0)) ;
    double strain_c[3] ;
    double strain_t[3] ;
    double strain_v = strain[0] + strain[4] + strain[8] ;
    double strain_tv = 0 ;
    double strain_cv = 0 ;
    double alpha_t = 0 ;
    double alpha_c = 0 ;
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      double stress_val = poisson / (1 - 2*poisson) * strain_v + strain_val[i] ;
      
      strain_t[i] = (stress_val > 0) ? stress_val : 0 ;
      strain_c[i] = (stress_val < 0) ? stress_val : 0 ;
      strain_tv += strain_t[i] ;
      strain_cv += strain_c[i] ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      strain_t[i] -= poisson / (1 + poisson) * strain_tv ;
      strain_c[i] -= poisson / (1 + poisson) * strain_cv ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      alpha_t += strain_t[i] * strain_pos[i] ;
      alpha_c += strain_c[i] * strain_pos[i] ;
    }
    
    alpha_t /= (strain_eq * strain_eq) ;
    alpha_c /= (strain_eq * strain_eq) ;
    
    {
      double dd = alpha_t * dd_t + alpha_c * dd_c ;
      double dd1 = (fabs(dd) > 0) ? 1./dd : 0 ;
    
      hm[0] = dd1 ;
    }
  }
  
  Damage_GetCriterionValue(damage) = crit ;
  return(crit) ;
}



double Damage_ReturnMappingMazars(Damage_t* damage,double* strain,double* d,double* hardv)
/** Mazars return mapping.
 *  Inputs are:
 *  the strains (strain), 
 *  the largest equivalent strain (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the Poisson's ratio (poisson),
 *  the limit elastic strain at the peak of uniaxial tensile stress test,
 *  the parameters of Mazars' model: A_t,B_t,A_c,B_c.
 * 
 *  On outputs, the following values are modified:
 *  the damage (d), 
 *  the largest equivalent strain reached in material history (kappa),
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double* cijkl  = Damage_GetStiffnessTensor(damage) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double kappa0  = Damage_GetStrainAtUniaxialTensileStrength(damage) ;
  double A_t     = Damage_GetA_t(damage) ;
  double B_t     = Damage_GetB_t(damage) ;
  double A_c     = Damage_GetA_c(damage) ;
  double B_c     = Damage_GetB_c(damage) ;

  double* strain_val ; //= Math_ComputePrincipalStresses(strain) ;
  double strain_dir[9] ;
  double strain_pos[3] ;
  double strain_eq ;
  double crit ;
  
  /* Principal values and principal directions */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      strain_dir[i] = strain[i] ;
    }
    
    strain_val = Math_ComputeRealEigenvaluesAndEigenvectorsOf3x3Matrix(strain_dir,'r') ;
  }
  
  /*
    Positive principal strains
  */
  {
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_pos[i] = Math_Max(strain_val[i],0) ;
    }
  }
  
  /*
    Equivalent tensile strain
  */
  {
    int    i ;
    
    strain_eq = 0 ;
    
    for(i = 0 ; i < 3 ; i++) {
      strain_eq += strain_pos[i]*strain_pos[i] ;
    }
    
    strain_eq = sqrt(strain_eq) ;
  }
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    
    crit = strain_eq - kappa ;
  }
  
  
  /*
    Return mapping: update damage and hardening variable
  */
  if(crit > 0.) {
    double kappa = strain_eq ;
    double d_t = 1 - kappa0 / kappa * (1 - A_t) - A_t * exp(-B_t * (kappa - kappa0)) ;
    double d_c = 1 - kappa0 / kappa * (1 - A_c) - A_c * exp(-B_c * (kappa - kappa0)) ;
    double strain_c[3] ;
    double strain_t[3] ;
    double strain_v = strain[0] + strain[4] + strain[8] ;
    double strain_tv = 0 ;
    double strain_cv = 0 ;
    double alpha_t = 0 ;
    double alpha_c = 0 ;
    int    i ;
    
    for(i = 0 ; i < 3 ; i++) {
      double stress_val = poisson / (1 - 2*poisson) * strain_v + strain_val[i] ;
      
      strain_t[i] = (stress_val > 0) ? stress_val : 0 ;
      strain_c[i] = (stress_val < 0) ? stress_val : 0 ;
      strain_tv += strain_t[i] ;
      strain_cv += strain_c[i] ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      strain_t[i] -= poisson / (1 + poisson) * strain_tv ;
      strain_c[i] -= poisson / (1 + poisson) * strain_cv ;
    }
    
    for(i = 0 ; i < 3 ; i++) {
      alpha_t += strain_t[i] * strain_pos[i] ;
      alpha_c += strain_c[i] * strain_pos[i] ;
    }
    
    alpha_t /= (strain_eq * strain_eq) ;
    alpha_c /= (strain_eq * strain_eq) ;


    d[0] = alpha_t * d_t + alpha_c * d_c ;
    hardv[0] = kappa ;
  }
  
  /* Current stiffness matrix */
  {
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    int i ;
    
    for(i = 0 ; i < 81 ; i++) {
      cijkl[i] = (1 - d[0]) * c[i] ;
    }
  }
  
  return(crit) ;
}







double Damage_ComputeFunctionGradientsMarigoJirasek(Damage_t* damage,const double* strain,const double* d,const double* hardv)
/** Energy release rate criterion
 *  Refs:
 *  [1] J.J. Marigo, Modelling of brittle and fatigue damage for elastic material 
 *  by growth of microvoids, Engineering Fracture Mechanics, 21(4):861-874, 1985.
 *  [2] M. Jirasek, B. Patzak, Consistent tangent stiffness for nonlocal damage models,
 *  Computers and Structures 80 (2002) 1279-1293.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest energy release rate (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the uniaxial tensile strength
 *  the fracture energy
 *  the width of the crack band
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  dgsde = derivative of the potential function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double* dfsde  = Damage_GetYieldFunctionGradient(damage) ;
  double* dgsde  = Damage_GetPotentialFunctionGradient(damage) ;
  double* hm     = Damage_GetHardeningModulus(damage) ;
  
  double ft      = Damage_GetUniaxialTensileStrength(damage) ;
  double Gf      = Damage_GetFractureEnergy(damage) ;
  double w       = Damage_GetCrackBandWidth(damage) ;
  double E       = Elasticity_GetYoungModulus(elasty) ;
  
  double eps0    = ft/E ;
  double gf      = Gf/w ;
  double epsf    = 0.5 * eps0 + gf/ft ;
  double kappa0  = 0.5 * E * eps0 * eps0 ;
  double kappaf  = 0.5 * E * epsf * epsf ;

  double crit ;
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    double energy = Elasticity_ComputeElasticEnergy(elasty,strain) ;
    /* G = energy release rate = - d(Free Energy)/d(d), i.e.
     * the opposite to the derivative of free energy wrt 
     * the damage variable at constant strains. 
     * Here G = 1/2 e:C0:e independant of the damage */
    double G = energy ;
    
    crit = G - kappa ;
  }
  
  /*
    Gradients
  */
  {
    int i ;
    
    Elasticity_ComputeStressTensor(elasty,strain,dfsde) ;
    
    for(i = 0 ; i < 9 ; i++) {
      dgsde[i] = dfsde[i] ;
    }
  }
  
  /* Hardening modulus: H = - d(crit)/d(d) = d(kappa)/d(d) */
  {
    double kappa = hardv[0] ;
    /* 
     * 1 - d = sqrt(kappa0/kappa) * exp(-(sqrt(kappa)-sqrt(kappa0)) / (sqrt(kappaf)-sqrt(kappa0)))
     */
    double mh = 0.5 * (1 - d[0]) / kappa * (1 + sqrt(kappa)/(sqrt(kappaf)-sqrt(kappa0))) ;
    
    hm[0] = 1/mh ;
  }
  
  Damage_GetCriterionValue(damage) = crit ;
  return(crit) ;
}





double Damage_ReturnMappingMarigoJirasek(Damage_t* damage,double* strain,double* d,double* hardv)
/** Energy release rate criterion
 *  Refs:
 *  [1] J.J. Marigo, Modelling of brittle and fatigue damage for elastic material 
 *  by growth of microvoids, Engineering Fracture Mechanics, 21(4):861-874, 1985.
 *  [2] M. Jirasek, B. Patzak, Consistent tangent stiffness for nonlocal damage models,
 *  Computers and Structures 80 (2002) 1279-1293.
 * 
 *  Inputs are:
 *  the strains (strain), 
 *  the largest energy release rate (kappa = hardv[0])
 * 
 *  Parameters are:
 *  the uniaxial tensile strength
 *  the fracture energy
 *  the width of the crack band
 * 
 *  On outputs, the following values are modified:
 *  dfsde = derivative of the yield function wrt strains
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Damage_GetElasticity(damage) ;
  double* cijkl  = Damage_GetStiffnessTensor(damage) ;
  
  double ft      = Damage_GetUniaxialTensileStrength(damage) ;
  double Gf      = Damage_GetFractureEnergy(damage) ;
  double w       = Damage_GetCrackBandWidth(damage) ;
  double E       = Elasticity_GetYoungModulus(elasty) ;
  
  double eps0    = ft/E ;
  double gf      = Gf/w ;
  double epsf    = 0.5 * eps0 + gf/ft ;
  double kappa0  = 0.5 * E * eps0 * eps0 ;
  double kappaf  = 0.5 * E * epsf * epsf ;

  double G ;
  double crit ;
    
  /*
    Criterion
  */
  {
    double kappa = hardv[0] ;
    double energy = Elasticity_ComputeElasticEnergy(elasty,strain) ;
    /* G = energy release rate, i.e. here
     * the opposite to the derivative of free energy wrt 
     * the damage variable at constant strains. */
    
    G = energy ;
    
    crit = G - kappa ;
  }
  
  
  /*
    Return mapping: update damage, hardening variable
   */
  if(crit > 0.) {
    double kappa = G ;
    double d1 = 1 - sqrt(kappa0 / kappa) * exp(-(sqrt(kappa)-sqrt(kappa0)) / (sqrt(kappaf)-sqrt(kappa0))) ;

    d[0] = d1 ;
    hardv[0] = kappa ;
  }
  
  /* Current stiffness matrix */
  {
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    int i ;
    
    for(i = 0 ; i < 81 ; i++) {
      cijkl[i] = (1 - d[0]) * c[i] ;
    }
  }
  
  return(crit) ;
}
