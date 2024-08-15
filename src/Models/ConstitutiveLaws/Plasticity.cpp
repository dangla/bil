#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Exception.h"
#include "Mry.h"
#include "Message.h"
#include "DataFile.h"
#include "Math_.h"
#include "Plasticity.h"
#include "BilExtraLibs.h"
#include "autodiff.h"


extern  Plasticity_SetModelProp_t Plasticity_ListOfSetModelProp ;

#define DEBUGGENERICRETURNMAPPING 0


static double* (Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm)(Plasticity_t*,const double*,const double*,const double*,double*) ;

static double* (Plasticity_ResidusOfGenericReturnMappingAlgorithm)(Plasticity_t*,const double*,const double*,const double*,const double*,const double*,double*) ;


Plasticity_t*  (Plasticity_Create)(void)
{
  Plasticity_t* plasty = (Plasticity_t*) Mry_New(Plasticity_t) ;
  
  
  /* Allocation of space for the code name of the model */
  {
    char* name = (char*) Mry_New(char[Plasticity_MaxLengthOfKeyWord]) ;
    
    Plasticity_GetCodeNameOfModel(plasty) = name ;
  }
  
  /* Allocation of space for the yield function gradient */
  {
    int n = (9 + Plasticity_MaxNbOfHardeningVariables)*Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n]) ;
    
    Plasticity_GetYieldFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the potential function gradient */
  {
    int n = (9 + Plasticity_MaxNbOfHardeningVariables)*Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n]) ;
    
    Plasticity_GetPotentialFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the hardening variable */
  {
    int n = 2*Plasticity_MaxNbOfHardeningVariables ;
    double* c = (double*) Mry_New(double[n]) ;
    
    Plasticity_GetHardeningVariable(plasty) = c ;
    Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) = c + Plasticity_MaxNbOfHardeningVariables ;
  }
  
  /* Allocation of space for the hardening modulus */
  {
    int n = Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n*n]) ;
    
    Plasticity_GetHardeningModulus(plasty) = c ;
  }
  
  /* Allocation of space for the criterion values */
  {
    double* c = (double*) Mry_New(double[Plasticity_MaxNbOfCriteria]) ;
    
    Plasticity_GetCriterionValue(plasty) = c ;
  }
  
  /* Allocation of space for the plastic multipliers */
  {
    double* c = (double*) Mry_New(double[Plasticity_MaxNbOfCriteria]) ;
    
    Plasticity_GetPlasticMultiplier(plasty) = c ;
  }
  
  /* Allocation of space for Fji*Cijkl */
  {
    int n = 9*Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n]) ;
    
    Plasticity_GetFjiCijkl(plasty) = c ;
  }
  
  /* Allocation of space for Cijkl*Glk */
  {
    int n = 9*Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n]) ;
    
    Plasticity_GetCijklGlk(plasty) = c ;
  }
  
  /* Allocation of space for Fji*Cijkl*Glk */
  {
    int n = Plasticity_MaxNbOfCriteria ;
    double* c = (double*) Mry_New(double[n*n]) ;
    
    Plasticity_GetFjiCijklGlk(plasty) = c ;
  }
  
  /* Allocation of space for the tangent stiffness tensor */
  {
    double* c = (double*) Mry_New(double[81]) ;
    
    Plasticity_GetTangentStiffnessTensor(plasty) = c ;
  }
  
  /* Elasticity */
  Plasticity_GetElasticity(plasty) = Elasticity_Create() ;
  
    
  /* The generic data */
  {
    Plasticity_GetGenericData(plasty) = NULL ;
  }
    
  /* The parameters (also part of the generic data) */
  {
    double* par = (double*) Mry_New(double[Plasticity_MaxNbOfParameters]) ;
    
    Plasticity_GetParameter(plasty) = par ;
      
    Plasticity_AppendData(plasty,Plasticity_MaxNbOfParameters,par,double,"Parameters") ;
  }

  /* Curves */
  {
    Plasticity_GetCurves(plasty) = Curves_Create(Plasticity_MaxNbOfCurves) ;
  }

  
  /* Space allocation for buffers */
  {
    Buffers_t* buf = Buffers_Create(Plasticity_SizeOfBuffer) ;
    
    Plasticity_GetBuffers(plasty) = buf ;
  }
  
  /* Default values */
  Plasticity_GetNbOfCriteria(plasty) = 1 ;
  Plasticity_GetNbOfHardeningVariables(plasty) = 0 ;
  Plasticity_GetNbOfNonHardeningVariables(plasty) = 0 ;
  
  return(plasty) ;
}



void  (Plasticity_Delete)(void* self)
{
  Plasticity_t* plasty = (Plasticity_t*) self ;
  
  {
    char* name = Plasticity_GetCodeNameOfModel(plasty) ;
    
    if(name) free(name) ;
  }
  
  {
    double* c = Plasticity_GetYieldFunctionGradient(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetPotentialFunctionGradient(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetHardeningVariable(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetCriterionValue(plasty) ;
    
    if(c) free(c) ;
      
  }
  
  {
    double* c = Plasticity_GetPlasticMultiplier(plasty) ;
    
    if(c) free(c) ;
      
  }
  
  {
    double* c = Plasticity_GetHardeningModulus(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetFjiCijkl(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetCijklGlk(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    
    if(c) free(c) ;
  }
  
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    if(elasty) {
      Elasticity_Delete(elasty) ;
      free(elasty) ;
      Plasticity_GetElasticity(plasty) = NULL ;
    }
  }
  
  {
    GenericData_t* gdat = Plasticity_GetGenericData(plasty) ;
    
    if(gdat) {
      GenericData_Delete(gdat) ;
      free(gdat) ;
      Plasticity_GetGenericData(plasty) = NULL ;
    }
  }
  
  {
    Curves_t* curves = Plasticity_GetCurves(plasty) ;
    
    if(curves) {
      /* The curves are deleted in Material_Delete */
      //Curves_Delete(curves) ;
      free(curves) ;
      Plasticity_GetCurves(plasty) = NULL ;
    }
  }
  
  {
    Buffers_t* buf = Plasticity_GetBuffers(plasty) ;

    if(buf) {
      Buffers_Delete(buf)  ;
      free(buf) ;
      Plasticity_GetBuffers(plasty) = NULL ;
    }
  }
}



void (Plasticity_Initialize)(Plasticity_t* plasty,const char* codename)
{
  int n_models = Plasticity_NbOfModels ;
  const char* modelnames[] = {Plasticity_ListOfNames} ;
  Plasticity_SetModelProp_t* xPlasticity_SetModelProp[] = {Plasticity_ListOfSetModelProp} ;
  int i = 0 ;
  
  while(i < n_models && strcmp(modelnames[i],codename)) i++ ;
    
  if(i < n_models) {
    Plasticity_CopyCodeName(plasty,modelnames[i]) ;
    Plasticity_GetSetModelProp(plasty) = xPlasticity_SetModelProp[i] ;
    Plasticity_SetModelProp(plasty) ;
    
    return ;
  } else {
    
    Message_Warning("No model named %s",codename) ;
  }

  return ;
}



#if 0
void (Plasticity_ScanProperties)(Plasticity_t* plasty,DataFile_t* datafile,int (*pm)(const char*))
/** Read the plastic properties in the stream file ficd */
{
  FILE* ficd = DataFile_GetFileStream(datafile) ;
  int    nd = Plasticity_GetNbOfProperties(plasty) ;
  short int    cont = 1 ;
  
  if(!ficd) return ;

  while(cont) {
    char   mot[Plasticity_MaxLengthOfKeyWord] ;
    char* line = DataFile_ReadLineFromCurrentFilePosition(datafile) ;
    
    if(!line) break ;
    
    sscanf(line," %[^= ]",mot) ;

    /* Reading some curves */
    if(!strncmp(mot,"Courbes",6) || !strncmp(mot,"Curves",5)) {
      Curves_t* curves = Plasticity_GetCurves(plasty) ;
      
      Curves_ReadCurves(curves,line) ;
      
      if(Curves_GetNbOfCurves(curves) > Plasticity_MaxNbOfCurves) {
        arret("Plasticity_ScanProperties (2) : trop de courbes") ;
      }
      
    /* Reading the material properties and storing through pm */
    } else if(pm) {
      char   *p = strchr(line,'=') ;

      /* We assume that this is a property as long as "=" is found */
      if(p) {
        int i = (*pm)(mot) ;
        
        if(i >= 0) {
        
          sscanf(p+1,"%lf",Plasticity_GetProperty(plasty) + i) ;
          nd = (nd > i + 1) ? nd : i + 1 ;
        
        } else {
        
          Message_RuntimeError("%s is not known",mot) ;
          
        }
      
      /* ... otherwise we stop reading */
      } else {
        
        /* go out */
        cont = 0 ;
      }
      
    } else {
      break ;
    }
    
  }

  Plasticity_GetNbOfProperties(plasty) = nd ;
  return ;
}
#endif





void (Plasticity_CopyTangentStiffnessTensor)(Plasticity_t* plasty,double* c)
/** Copy the 4th rank tangent stiffness tensor in c. */
{
  double* ct = Plasticity_GetTangentStiffnessTensor(plasty) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = ct[i] ;
    }
  }
}



int (Plasticity_UpdateElastoplasticTensor)(Plasticity_t* plasty,double* c)
/** Update the 4th rank elastoplastic tensor in c.
 *  On input c should point to the elastic stiffness tensor. 
 *  On output c is updated to the tangent elastoplastic stiffness tensor.
 *  Other inputs, included in plasty, are: 
 *  dfsds is the yield function gradient
 *  dgsds is the potential function gradient
 *  hm    is the hardening modulus
 *  Other outputs, saved in plasty, are:
 *  fc(k,l) = dfsds(j,i) * C^el(i,j,k,l)
 *  cg(i,j) = C^el(i,j,k,l) * dgsds(l,k)
 *  fcg     = dfsds(j,i) * C^el(i,j,k,l) * dgsds(l,k)
 *  Tensor c is then updated as 
 *  C(i,j,k,l) = C^el(i,j,k,l) + cg(i,j) * fc(k,l) / det
 *  with det = hm + dfsds(j,i) * C^el(i,j,k,l) * dgsds(l,k)
 *  Return the inverse of det: 1/det */
{
#define CEP(i,j)  ((c)[(i)*9+(j)])
#define CEL(i,j)  ((c)[(i)*9+(j)])
//#define CEL(i,j)  ((cel)[(i)*9+(j)])
  //Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  //double* cel    = Elasticity_GetStiffnessTensor(elasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double* fc     = Plasticity_GetFjiCijkl(plasty) ;
  double* cg     = Plasticity_GetCijklGlk(plasty) ;
  double  fcg    = 0 ;
  double  det ;
  
  if(Plasticity_GetNbOfCriteria(plasty) > 1) {
    Message_FatalError("Plasticity_UpdateElastoplasticTensor: not available yet!") ;
  }
  
  /* Tangent elastoplastic stiffness tensor */
  {
      int i ;
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        fc[i] = 0. ;
        cg[i] = 0. ;
          
        for(j = 0 ; j < 9 ; j++) {
              
          fc[i] += dfsds[j]*CEL(j,i) ;
          cg[i] += CEL(i,j)*dgsds[j] ;
          fcg   += dfsds[i]*CEL(i,j)*dgsds[j] ;
        }
      }
        
      Plasticity_GetFjiCijklGlk(plasty)[0] = fcg ;
      
      det = hm[0] + fcg ;
        
      if(det > 0.) {
        det = 1./det ;
      } else {
            
        Message_Direct("\n") ;
        Message_Direct("Plasticity_UpdateElastoplasticTensor:") ;
        Message_Direct("\n") ;
        Message_Direct("Something wrong happened\n") ;
        
        printf("\n") ;
        printf("dF = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dfsds[i]) ;
        }
        printf("\n") ;
        printf("dG = ") ;
        for(i = 0 ; i < 9 ; i++) {
          printf(" %e",dgsds[i]) ;
        }
        printf("\n") ;
        printf("hm = %e\n",hm[0]) ;
        printf("dF(j,i) * C(i,j,k,l) * dG(l,k) = %e\n",fcg) ;
        printf("hm + dF(j,i) * C(i,j,k,l) * dG(l,k) = %e\n",det) ;
        printf("\n") ;
        
        //Exception_BackupAndTerminate ;
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          CEP(i,j) = CEL(i,j) - cg[i]*fc[j]*det ;
        }
      }
  }
  
  return(0) ;
#undef CEP
#undef CEL
}




void (Plasticity_PrintTangentStiffnessTensor)(Plasticity_t* plasty)
/** Print the 4th rank tangent elastoplastic tensor.
 **/
{
  double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
  
  printf("\n") ;
  printf("4th rank elastoplastic tensor:\n") ;
  
  Math_PrintStiffnessTensor(c) ;
}


#if 1
template<>
double* (Plasticity_DerivativeOfYieldFunction<Plasticity_YieldFunction_t*>)(Plasticity_YieldFunction_t* yieldfunction,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the yield functions w.r.t. the stresses 
 *  and the hardening variables. 
 *  The output is an array of (9+nhardv)*ncrit components. 
 *  For each criterium the 9 first components are the derivatives of the 
 *  yield function w.r.t. the stresses and the following components are 
 *  the derivative w.r.t. the hardening variables.
 */
{
  double dstress = Plasticity_GetTypicalSmallIncrementOfStress(plasty) ;
  double* dhardv = Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  size_t SizeNeeded = nflows*ncrit*(sizeof(double)) ;
  double* dyield = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  if(dstress == 0) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  /* k ranges from 0 to ncrit-1 (index of the yield function) */
  #define DYIELD(k)          (dyield + nflows*(k))
  #define DYIELD_STRESS(k)   (DYIELD(k))
  #define DYIELD_HARDV(k)    (DYIELD(k) + 9)
  {
    double* yield0 = yieldfunction(plasty,stress,hardv) ;
    double stress1[9] ;
    double hardv1[Plasticity_MaxNbOfHardeningVariables] ;
    int i ;
    
    for(i = 0 ; i < nthardv ; i++) {
      hardv1[i] = hardv[i] ;
    }
    
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] = stress[i] ;
    }
    
    /* Derivatives with respect to stresses */
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] += dstress ;
      
      {
        double* yield1 = yieldfunction(plasty,stress1,hardv) ;
        int k ;
    
        for(k = 0 ; k < ncrit ; k++) {
          double* dyield_stress = DYIELD_STRESS(k) ;
          
          dyield_stress[i] = (yield1[k] - yield0[k])/dstress ;
        }
        
        Plasticity_FreeBufferFrom(plasty,yield1) ;
      }
      
      stress1[i] = stress[i] ;
    }
    
    /* Derivatives with respect to hardening variables */
    for(i = 0 ; i < nhardv ; i++) {
      
      if(dhardv[i] == 0) {
        arret("Plasticity_DerivativeOfYieldFunction") ;
      }
      
      hardv1[i] += dhardv[i] ;
      
      {
        double* yield1 = yieldfunction(plasty,stress,hardv1) ;
        int k ;
    
        for(k = 0 ; k < ncrit ; k++) {
          double* dyield_hardv = DYIELD_HARDV(k) ;
          
          dyield_hardv[i] = (yield1[k] - yield0[k])/dhardv[i] ;
        }
        
        Plasticity_FreeBufferFrom(plasty,yield1) ;
      }
      
      hardv1[i] = hardv[i] ;
    }
        
    Plasticity_FreeBufferFrom(plasty,yield0) ;
  }
  #undef DYIELD_STRESS
  #undef DYIELD_HARDV
  #undef DYIELD
  
  return(dyield) ;
}
#endif


#ifdef HAVE_AUTODIFF
template<>
double* (Plasticity_DerivativeOfYieldFunction<Plasticity_YieldFunctionDual_t*>)(Plasticity_YieldFunctionDual_t* yieldfunction,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the yield functions w.r.t. the stresses 
 *  and the hardening variables using autodiff.
 *  The output is an array of (9+nhardv)*ncrit components. 
 *  For each criterium the 9 first components are the derivatives of the 
 *  yield function w.r.t. the stresses and the following components are 
 *  the derivative w.r.t. the hardening variables.
 */
{
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  size_t SizeNeeded = nflows*ncrit*(sizeof(double)) ;
  double* dyield = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  /* k ranges from 0 to ncrit-1 (index of the yield function) */
  #define DYIELD(k)          (dyield + nflows*(k))
  #define DYIELD_STRESS(k)   (DYIELD(k))
  #define DYIELD_HARDV(k)    (DYIELD(k) + 9)
  {
    Eigen::MatrixXd yieldg;
    
    {
      std::vector<real> yield(ncrit);
      std::vector<real> stress1(9) ;
      std::vector<real> hardv1(nthardv) ;
      int i ;
    
      for(i = 0 ; i < nthardv ; i++) {
        hardv1[i] = hardv[i] ;
      }
    
      for(i = 0; i < 9; i++) {
        stress1[i] = stress[i] ;
      }
 
      /* We use a lambda function */
      auto yf1 = [](real* (*yf)(Plasticity_t*,const real*,const real*),Plasticity_t* p,const std::vector<real> s,const std::vector<real> h)->std::vector<real>{
        real* y = yf(p,&(s[0]),&(h[0])) ;
        int n  = Plasticity_GetNbOfCriteria(p) ;
        std::vector<real> y1(n) ;
        
        for(int i = 0 ; i < n ; i++) {
          y1[i] = y[i] ;
        }
        
        Plasticity_FreeBufferFrom(p,y) ;
        
        return(y1) ;
      };
    
      yieldg = jacobian(yf1, wrt(stress1,hardv1), at(yieldfunction,plasty,stress1,hardv1),yield);
    }
    
    /* Derivatives with respect to stresses */
    {
      int k ;
    
      for(k = 0 ; k < ncrit ; k++) {
        double* dyield_stress = DYIELD_STRESS(k) ;
        int i ;
    
        for(i = 0 ; i < 9 ; i++) {
          dyield_stress[i] = yieldg(k,i) ;
        }
      }
    }
    
    /* Derivatives with respect to hardening variables */
    {
      int k ;
    
      for(k = 0 ; k < ncrit ; k++) {
        double* dyield_hardv = DYIELD_HARDV(k) ;
        int i ;
        
        for(i = 0 ; i < nhardv ; i++) {
          dyield_hardv[i] = yieldg(k,9+i) ;
        }
      }
    }
  }
  #undef DYIELD_STRESS
  #undef DYIELD_HARDV
  #undef DYIELD
  
  return(dyield) ;
}
#endif


#if 1
template<>
double* (Plasticity_DerivativeOfYieldFunction<Plasticity_YieldFunction_ftor>)(Plasticity_YieldFunction_ftor yieldfunction,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the yield functions w.r.t. the stresses 
 *  and the hardening variables using autodiff or not.
 *  The output is an array of (9+nhardv)*ncrit components. 
 *  For each criterium the 9 first components are the derivatives of the 
 *  yield function w.r.t. the stresses and the following components are 
 *  the derivative w.r.t. the hardening variables.
 */
{
  Plasticity_YieldFunction_t*      yf  = yieldfunction.YieldFunction;
  Plasticity_YieldFunctionDual_t*  yfd = yieldfunction.YieldFunctionDual;
  
  if(0) {
  #ifdef HAVE_AUTODIFF
  } else if(yfd) {
    return(Plasticity_DerivativeOfYieldFunction(yfd,plasty,stress,hardv));
  #endif
  } else if(yf) {
    return(Plasticity_DerivativeOfYieldFunction(yf,plasty,stress,hardv));
  } else {
    arret("Plasticity_DerivativeOfYieldFunction:") ;
  }
}
#endif


#if 1
template<>
double* (Plasticity_DerivativeOfFlowRules<Plasticity_FlowRules_t*>)(Plasticity_FlowRules_t* flowrules,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the flow rules w.r.t. the stresses 
 *  and the hardening variables.
 *  The outputs are composed of ncrit matrices (9+nhardv)x(9+nhardv). 
 *  The components of each matrix are stored in the following order:
 * 
 *  K0 to Kn then A0 to An etc.. with n = nhardv:
 * 
 *  | K0(9x9) K1(9x1) K2(9x1) ... | 
 *  | A0(1x9) A1(1x1) A2(1x1) ... | 
 *  | B0(1x9) B1(1x1) B2(1x1) ... | 
 *  | ........................... |
 * 
 *   K0    = derivative of the plastic strain flow w.r.t. the stresses
 *   Ki    = derivative of the plastic strain flow w.r.t. the ith hardening variable
 *   A0    = derivative of the first hardening flow w.r.t. the stresses
 *   Ai    = derivative of the first hardening flow w.r.t. the ith hardening variables
 *   B0    = derivative of the second hardening flow w.r.t. the stresses
 *   Bi    = derivative of the second hardening flow w.r.t. the ith hardening variables
 *   etc...
 */
{
  double dstress = Plasticity_GetTypicalSmallIncrementOfStress(plasty) ;
  double* dhardv = Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  int nflows2 = nflows*nflows ;
  size_t SizeNeeded = nflows2*ncrit*(sizeof(double)) ;
  double* dflow = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
  
  if(dstress == 0) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
    
  /* k ranges from 0 to ncrit-1 (index of the yield mechanism) */
  #define DFLOW(k)                (dflow + nflows2*(k))
  #define DFLOWSTRAIN_STRESS(k)   (DFLOW(k))
  #define DFLOWSTRAIN_HARDV(k)    (DFLOW(k) + 81)
  #define DFLOWHARDV_STRESS(k)    (DFLOW(k) + 9*nflows)
  #define DFLOWHARDV_HARDV(k)     (DFLOW(k) + 9*nflows + 9)
  {
    double* flow0 = flowrules(plasty,stress,hardv) ;
    double stress1[9] ;
    double hardv1[Plasticity_MaxNbOfHardeningVariables] ;
    int i ;
    
    for(i = 0 ; i < nthardv ; i++) {
      hardv1[i] = hardv[i] ;
    }
    
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] = stress[i] ;
    }
    
    /* Derivatives with respect to stresses */
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] += dstress ;
      
      {
        double* flow1 = flowrules(plasty,stress1,hardv) ;
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_stress =  DFLOWSTRAIN_STRESS(k) ;
          double* dflowhardv_stress =  DFLOWHARDV_STRESS(k) ;
          double* flowk0 = flow0 + nflows*k ;
          double* flowk1 = flow1 + nflows*k ;
          int j ;
    
          for(j = 0 ; j < 9 ; j++) {
            dflowstrain_stress[j*9+i] = (flowk1[j] - flowk0[j])/dstress ;
          }
    
          for(j = 0 ; j < nhardv ; j++) {
            dflowhardv_stress[j*nflows+i] = (flowk1[9+j] - flowk0[9+j])/dstress ;
          }
        }
        
        Plasticity_FreeBufferFrom(plasty,flow1) ;
      }
      
      stress1[i] = stress[i] ;
    }
    
    /* Derivatives with respect to hardening variables */
    for(i = 0 ; i < nhardv ; i++) {
      
      if(dhardv[i] == 0) {
        arret("Plasticity_DerivativeOfFlowRules") ;
      }
      
      hardv1[i] += dhardv[i] ;
      
      {
        double* flow1 = flowrules(plasty,stress,hardv1) ;
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_hardv = DFLOWSTRAIN_HARDV(k) ;
          double* dflowhardv_hardv = DFLOWHARDV_HARDV(k) ;
          double* flowk0 = flow0 + nflows*k ;
          double* flowk1 = flow1 + nflows*k ;
          int j ;
    
          for(j = 0 ; j < 9 ; j++) {
            dflowstrain_hardv[i*9+j] = (flowk1[j] - flowk0[j])/dhardv[i] ;
          }
    
          for(j = 0 ; j < nhardv ; j++) {
            dflowhardv_hardv[j*nflows+i] = (flowk1[9+j] - flowk0[9+j])/dhardv[i] ;
          }
        }
        
        Plasticity_FreeBufferFrom(plasty,flow1) ;
      }
      
      hardv1[i] = hardv[i] ;
    }
    
    Plasticity_FreeBufferFrom(plasty,flow0) ;
  }
  #undef DFLOWSTRAIN_STRESS
  #undef DFLOWSTRAIN_HARDV
  #undef DFLOWHARDV_STRESS
  #undef DFLOWHARDV_HARDV
  #undef DFLOW
  
  return(dflow) ;
}
#endif


#ifdef HAVE_AUTODIFF
template<>
double* (Plasticity_DerivativeOfFlowRules<Plasticity_FlowRulesDual_t*>)(Plasticity_FlowRulesDual_t* flowrules,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the flow rules w.r.t. the stresses 
 *  and the hardening variables using autodiff.
 *  The outputs are composed of ncrit matrices (9+nhardv)x(9+nhardv). 
 *  The components of each matrix are stored in the following order:
 * 
 *  K0 to Kn then A0 to An etc.. with n = nhardv:
 * 
 *  | K0(9x9) K1(9x1) K2(9x1) ... | 
 *  | A0(1x9) A1(1x1) A2(1x1) ... | 
 *  | B0(1x9) B1(1x1) B2(1x1) ... | 
 *  | ........................... |
 * 
 *   K0    = derivative of the plastic strain flow w.r.t. the stresses
 *   Ki    = derivative of the plastic strain flow w.r.t. the ith hardening variable
 *   A0    = derivative of the first hardening flow w.r.t. the stresses
 *   Ai    = derivative of the first hardening flow w.r.t. the ith hardening variables
 *   B0    = derivative of the second hardening flow w.r.t. the stresses
 *   Bi    = derivative of the second hardening flow w.r.t. the ith hardening variables
 *   etc...
 */
{
  double dstress = Plasticity_GetTypicalSmallIncrementOfStress(plasty) ;
  double* dhardv = Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  int nflows2 = nflows*nflows ;
  size_t SizeNeeded = nflows2*ncrit*(sizeof(double)) ;
  double* dflow = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
    
  /* k ranges from 0 to ncrit-1 (index of the yield mechanism) */
  #define DFLOW(k)                (dflow + nflows2*(k))
  #define DFLOWSTRAIN_STRESS(k)   (DFLOW(k))
  #define DFLOWSTRAIN_HARDV(k)    (DFLOW(k) + 81)
  #define DFLOWHARDV_STRESS(k)    (DFLOW(k) + 9*nflows)
  #define DFLOWHARDV_HARDV(k)     (DFLOW(k) + 9*nflows + 9)
  {
    Eigen::MatrixXd flowg;
    int ntflows = 9 + nthardv ;
    int i ;
    
    {
      std::vector<real> flow(9+nthardv);
      std::vector<real> stress1(9) ;
      std::vector<real> hardv1(nthardv) ;
    
      for(i = 0 ; i < nthardv ; i++) {
        hardv1[i] = hardv[i] ;
      }
    
      for(i = 0; i < 9; i++) {
        stress1[i] = stress[i] ;
      }
 
      auto fr1 = [](real* (*fr)(Plasticity_t*,const real*,const real*),Plasticity_t* p,const std::vector<real> s,const std::vector<real> h) {
        real* f = fr(p,&(s[0]),&(h[0])) ;
        int nhardv = Plasticity_GetNbOfHardeningVariables(p) ;
        int nflows = 9 + nhardv ;
        int ncrit  = Plasticity_GetNbOfCriteria(p) ;
        int n = nflows*ncrit ;
        std::vector<real> f1(n) ;
        
        for(int i = 0 ; i < n ; i++) {
          f1[i] = f[i] ;
        }
        
        Plasticity_FreeBufferFrom(p,f) ;
        
        return(f1) ;
      };
    
      flowg = jacobian(fr1, wrt(stress1,hardv1), at(flowrules,plasty,stress1,hardv1),flow);
    }
    
    /* Derivatives with respect to stresses */
    for(i = 0 ; i < 9 ; i++) {
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_stress = DFLOWSTRAIN_STRESS(k) ;
          double* dflowhardv_stress  = DFLOWHARDV_STRESS(k) ;
          int j ;
    
          for(j = 0 ; j < 9 ; j++) {
            dflowstrain_stress[j*9+i] = flowg(j+ntflows*k,i) ;
          }
    
          for(j = 0 ; j < nhardv ; j++) {
            dflowhardv_stress[j*nflows+i] = flowg(j+9+ntflows*k,i) ;
          }
        }
      }
    }
    
    /* Derivatives with respect to hardening variables */
    for(i = 0 ; i < nhardv ; i++) {
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_hardv = DFLOWSTRAIN_HARDV(k) ;
          double* dflowhardv_hardv  = DFLOWHARDV_HARDV(k) ;
          int j ;
    
          for(j = 0 ; j < 9 ; j++) {
            dflowstrain_hardv[i*9+j] = flowg(j+ntflows*k,9+i) ;
          }
    
          for(j = 0 ; j < nhardv ; j++) {
            dflowhardv_hardv[j*nflows+i] = flowg(j+9+ntflows*k,9+i) ;
          }
        }
      }
    }
  }
  #undef DFLOWSTRAIN_STRESS
  #undef DFLOWSTRAIN_HARDV
  #undef DFLOWHARDV_STRESS
  #undef DFLOWHARDV_HARDV
  #undef DFLOW
  
  return(dflow) ;
}
#endif


#if 1
template<>
double* (Plasticity_DerivativeOfFlowRules<Plasticity_FlowRules_ftor>)(Plasticity_FlowRules_ftor flowrules,Plasticity_t* plasty,const double* stress,const double* hardv)
/** Return the derivative of the flow rules w.r.t. the stresses 
 *  and the hardening variables using autodiff or not.
 *  The outputs are composed of ncrit matrices (9+nhardv)x(9+nhardv). 
 *  The components of each matrix are stored in the following order:
 * 
 *  K0 to Kn then A0 to An etc.. with n = nhardv:
 * 
 *  | K0(9x9) K1(9x1) K2(9x1) ... | 
 *  | A0(1x9) A1(1x1) A2(1x1) ... | 
 *  | B0(1x9) B1(1x1) B2(1x1) ... | 
 *  | ........................... |
 * 
 *   K0    = derivative of the plastic strain flow w.r.t. the stresses
 *   Ki    = derivative of the plastic strain flow w.r.t. the ith hardening variable
 *   A0    = derivative of the first hardening flow w.r.t. the stresses
 *   Ai    = derivative of the first hardening flow w.r.t. the ith hardening variables
 *   B0    = derivative of the second hardening flow w.r.t. the stresses
 *   Bi    = derivative of the second hardening flow w.r.t. the ith hardening variables
 *   etc...
 */
{
  Plasticity_FlowRules_t*     fr = flowrules.FlowRules;
  Plasticity_FlowRulesDual_t* frd = flowrules.FlowRulesDual;
  
  if(0) {
  #ifdef HAVE_AUTODIFF
  } else if(frd) {
    return(Plasticity_DerivativeOfFlowRules(frd,plasty,stress,hardv));
  #endif
  } else if(fr) {
    return(Plasticity_DerivativeOfFlowRules(fr,plasty,stress,hardv));
  } else {
    arret("Plasticity_DerivativeOfFlowRules:") ;
  }
}
#endif


#if 1
double* (Plasticity_ResidusOfGenericReturnMappingAlgorithm)(Plasticity_t* plasty,const double* stress,const double* hardv,const double* stress_t,const double* hardv_n,const double* zeta,double* residu)
/** Compute the residu vector used for the return mapping algorithm. 
 *  
 *  On input residu should point to an array of (9+nhardv+ncrit) doubles.
 * 
 *  On output the (6+nhardv+ncrit) first components of the residu vector 
 *  are associated to
 *    - the plastic strain flow rule
 *    - the hardening flow rule
 *    - the yield functions
 * 
 *  Return the pointer to the residu
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  
  if(!yieldfunction) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  if(!flowrules) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  
  {
    double* invc = Elasticity_GetComplianceTensor(elasty) ;
      
    /* k ranges from 0 to ncrit-1 (index of the yield mechanism) */
    #define FLOW(k)                   (flow + nflows*(k))
    #define FLOWSTRAIN(k)             (FLOW(k))
    #define FLOWHARDV(k)              (FLOW(k) + 9)
    #define RESIDUSTRAIN              (residu)
    #define RESIDUHARDV               (residu+6)
    #define RESIDUYIELD               (residu+6+nhardv)
    {
      double* yield = yieldfunction(plasty,stress,hardv) ;
      double* flow = flowrules(plasty,stress,hardv) ;
      
      /* Equation of the plastic strain rules */
      {
        double dstress[9] ;
        int i ;
        
        for(i = 0 ; i < 9 ; i++) {
          dstress[i] = (stress[i] - stress_t[i]) ;
        }
        
        for(i = 0 ; i < 9 ; i++) {
          int j ;
          
          RESIDUSTRAIN[i] = 0 ;
          
          for(j = 0 ; j < 9 ; j++) {
            RESIDUSTRAIN[i] += invc[9*i + j] * dstress[j] ;
          }
        }
      }
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* flowstrain = FLOWSTRAIN(k) ;
          double lambda = Math_Max(zeta[k],0) ;
          
          if(lambda > 0) {
            int i ;
            
            for(i = 0 ; i < 9 ; i++) {
              RESIDUSTRAIN[i] += lambda * flowstrain[i] ;
            }
          }
        }
        
        Elasticity_ConvertStressTensorInto6TermStressVector(RESIDUSTRAIN) ;
      }
      
      /* Equation of the hardening flow rules */
      {
        {
          int i ;
        
          for(i = 0 ; i < nhardv ; i++) {
            RESIDUHARDV[i] = hardv[i] - hardv_n[i] ;
          }
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* flowhardv = FLOWHARDV(k) ;
          double lambda = Math_Max(zeta[k],0) ;
          
          if(lambda > 0) {
            int i ;
        
            for(i = 0 ; i < nhardv ; i++) {
              RESIDUHARDV[i] -= lambda * flowhardv[i] ;
            }
          }
        }
      }
      
      /* Equation of the yield criteria */
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double yval = Math_Min(zeta[k],0) ;
        
          RESIDUYIELD[k] = yield[k] - yval ;
        }
      }
      
      /* Change the sign */
      {
        int nmat = 6 + nhardv + ncrit ;
        int i ;
        
        for(i = 0 ; i < nmat ; i++) {
          residu[i] = - residu[i] ;
        }
      }
      
      Plasticity_FreeBufferFrom(plasty,flow) ;
    }
    #undef FLOWSTRAIN
    #undef FLOWHARDV
    #undef RESIDUSTRAIN
    #undef RESIDUHARDV
    #undef RESIDUYIELD
  }
  
  return(residu) ;
}
#endif



#if 1
double* (Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm)(Plasticity_t* plasty,const double* stress,const double* hardv,const double* zeta,double* matrix)
/** Compute the jacobian matrix of the function residus used for 
 *  the return mapping algorithm.
 * 
 *  On input matrix should point to an array of 
 *  (6+nhardv+ncrit)x(6+nhardv+ncrit) doubles.
 * 
 *  On output the (6+nhardv+ncrit)x(6+nhardv+ncrit) matrix
 *  is associated to
 *    - the 6 stresses
 *    - the nhardv hardening variables
 *    - the ncrit zeta variables (plastic multipliers or yield functions)
 * 
 *  Return the pointer to the matrix
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  //Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  //Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  //Plasticity_YieldFunctionDual_t* yieldfunction = Plasticity_GetYieldFunctionDual(plasty) ;
  Plasticity_YieldFunction_ftor yieldfunction = Plasticity_GetYieldFunctionFtor(plasty) ;
  //Plasticity_FlowRulesDual_t* flowrules = Plasticity_GetFlowRulesDual(plasty) ;
  Plasticity_FlowRules_ftor flowrules = Plasticity_GetFlowRulesFtor(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  int nflows2 = nflows*nflows ;
  
  #if 0
  if(!yieldfunction) {
    arret("Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm") ;
  }
  #endif
  
  #if 0
  if(!flowrules) {
    arret("Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm") ;
  }
  #endif
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm") ;
  }
  
        
  /* k ranges from 0 to ncrit-1 (index of the yield mechanism) */
  {
    #define DYIELD(k)                 (dyield + nflows*(k))
    #define FLOW(k)                   (flow + nflows*(k))
    #define DFLOW(k)                  (dflow + nflows2*(k))
    #define DYIELD_STRESS(k)          (DYIELD(k))
    #define DYIELD_HARDV(k)           (DYIELD(k) + 9)
    #define FLOWSTRAIN(k)             (FLOW(k))
    #define FLOWHARDV(k)              (FLOW(k) + 9)
    #define DFLOWSTRAIN_STRESS(k)     (DFLOW(k))
    #define DFLOWSTRAIN_HARDV(k)      (DFLOW(k) + 81)
    #define DFLOWHARDV_STRESS(k)      (DFLOW(k) + 9*nflows)
    #define DFLOWHARDV_HARDV(k)       (DFLOW(k) + 9*nflows + 9)
    #define MATRIX(i,j)               matrix[(i)*nmat+(j)]
    {
      double* dyield = Plasticity_DerivativeOfYieldFunction(yieldfunction,plasty,stress,hardv) ;
      Plasticity_FlowRules_t* flowrules_double = Plasticity_GetFlowRules(plasty) ;
      double* flow   = flowrules_double(plasty,stress,hardv) ;
      double* dflow  = Plasticity_DerivativeOfFlowRules(flowrules,plasty,stress,hardv) ;
      int nmat = 6 + nhardv + ncrit ;
      
      {
        int i ;
        
        for(i = 0 ; i < nmat*nmat ; i++) {
          matrix[i] = 0 ;
        }
      }
      
      {
        double c1[81] ;
        double* invc = Elasticity_CopyComplianceTensor(elasty,c1) ;
        double* invc66 = Elasticity_ConvertStiffnessMatrixInto6x6Matrix(invc) ;
        int i ;
        
        for(i = 0 ; i < 6 ; i++) {
          int j ;
          
          for(j = 0 ; j < 6 ; j++) {
            MATRIX(i,j) = invc66[6*i+j] ;
          }
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_stress = DFLOWSTRAIN_STRESS(k) ;
          double lambda = Math_Max(zeta[k],0) ;
          
          if(lambda > 0) {
            int i ;
      
            Elasticity_ConvertStiffnessMatrixInto6x6Matrix(dflowstrain_stress) ;
        
            for(i = 0 ; i < 6 ; i++) {
              int j ;
          
              for(j = 0 ; j < 6 ; j++) {
                MATRIX(i,j) += lambda * dflowstrain_stress[6*i+j] ;
              }
            }
          }
        }
      }
      
      
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowstrain_hardv = DFLOWSTRAIN_HARDV(k) ;
          double lambda = Math_Max(zeta[k],0) ;
        
          if(lambda > 0) {
            int i ;
            
            for(i = 0 ; i < nhardv ; i++) {
              int j ;
              double* dflowk = dflowstrain_hardv + 9*i ;
          
              Elasticity_ConvertStressTensorInto6TermStressVector(dflowk) ;
          
              for(j = 0 ; j < 6 ; j++) {
                MATRIX(j,6+i) += lambda * dflowk[j] ;
              }
            }
          }
        }
      }
      
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* flowstrain = FLOWSTRAIN(k) ;
          double dlambda = (zeta[k] > 0) ? 1 : 0 ;
          
          if(dlambda > 0) {
            int i ;
            
            Elasticity_ConvertStressTensorInto6TermStressVector(flowstrain) ;
        
            for(i = 0 ; i < 6 ; i++) {
              MATRIX(i,6+nhardv+k) = dlambda * flowstrain[i] ;
            }
          }
        }
      }
      
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowhardv_stress = DFLOWHARDV_STRESS(k) ;
          double lambda = Math_Max(zeta[k],0) ;
          
          if(lambda > 0) {
            int i ;
        
            for(i = 0 ; i < nhardv ; i++) {
              int j ;
              double* dflowk = dflowhardv_stress + i*nflows ;
          
              Elasticity_ConvertStressTensorInto6TermStressVector(dflowk) ;
          
              for(j = 0 ; j < 6 ; j++) {
                MATRIX(6+i,j) = - lambda * dflowk[j] ;
              }
            }
          }
        }
      }

      {
        int i ;
        
        for(i = 0 ; i < nhardv ; i++) {
          MATRIX(6+i,6+i) = 1 ;
        }
      }

      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dflowhardv_hardv = DFLOWHARDV_HARDV(k) ;
          double lambda = Math_Max(zeta[k],0) ;
          
          if(lambda > 0) {
            int i ;
        
            for(i = 0 ; i < nhardv ; i++) {
              double* dflowk = dflowhardv_hardv + i*nflows ;
              int j ;
          
              for(j = 0 ; j < nhardv ; j++) {
                MATRIX(6+i,6+j) -= lambda * dflowk[j] ;
              }
            }
          }
        }
      }
      
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* flowhardv = FLOWHARDV(k) ;
          double dlambda = (zeta[k] > 0) ? 1 : 0 ;
          
          if(dlambda > 0) {
            int i ;
        
            for(i = 0 ; i < nhardv ; i++) {
              MATRIX(6+i,6+nhardv+k) = - dlambda * flowhardv[i] ;
            }
          }
        }
      }
      
      
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dyield_stress = DYIELD_STRESS(k) ;
          int j ;
        
          Elasticity_ConvertStressTensorInto6TermStressVector(dyield_stress) ;
        
          for(j = 0 ; j < 6 ; j++) {
            MATRIX(6+nhardv+k,j) = dyield_stress[j] ;
          }
        }
      }
      
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double* dyield_hardv = DYIELD_HARDV(k) ;
          int i ;
        
          for(i = 0 ; i < nhardv ; i++) {
            MATRIX(6+nhardv+k,6+i) = dyield_hardv[i] ;
          }
        }
      }
      
        
      {
        int k ;
        
        for(k = 0 ; k < ncrit ; k++) {
          double dyval = (zeta[k] > 0) ? 0 : 1 ;
        
          MATRIX(6+nhardv+k,6+nhardv+k) = - dyval ;
        }
      }
      
      Plasticity_FreeBufferFrom(plasty,dyield) ;
    }
    #undef DYIELD
    #undef FLOW
    #undef DFLOW
    #undef FLOWSTRAIN
    #undef FLOWHARDV
    #undef DFLOWSTRAIN_STRESS
    #undef DFLOWSTRAIN_HARDV
    #undef DFLOWHARDV_STRESS
    #undef DFLOWHARDV_HARDV
    #undef DYIELD_HARDV
    #undef DYIELD_STRESS
    #undef MATRIX
  }
  
  return(matrix) ;
}
#endif



#if 1
double* (Plasticity_GenericReturnMapping)(Plasticity_t* plasty,double* stress,double* strain_p,double* hardv)
/** Return mapping algorithm in a generic way
 * 
 *  On inputs: 
 *    - "stress" points to the current elastic trial stress tensor
 *    - "hardv" points to the hardening variables obtained at the
 *      previous time step
 *    - "strain_p" points to the plastic strains obtained at the
 *      previous time step
 * 
 *  On outputs, the following values are modified:
 *    - the stresses ("stress"), 
 *    - the plastic strains ("strain_p"), 
 *    - the hardening variables ("hardv").
 * 
 *  Return the value of the yield functions. 
 **/
{
  double* crit  = Plasticity_GetCriterionValue(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nflows = 9 + nhardv ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericReturnMapping") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_GenericReturnMapping") ;
  }
  
  {
    double stress_t[9] ;
    double hardv_n[Plasticity_MaxNbOfHardeningVariables] ;
    double zeta[Plasticity_MaxNbOfCriteria] ;
    int iter = 0 ;
    int niter = 20 ;
    int convergencenotattained = 1 ;
  
    {
      int i ;
    
      for(i = 0 ; i < 9 ; i++) {
        stress_t[i] = stress[i] ;
      }
    
      for(i = 0 ; i < nthardv ; i++) {
        hardv_n[i] = hardv[i] ;
      }
    
      for(i = 0 ; i < ncrit ; i++) {
        zeta[i] = 0 ;
      }
    }
    
    while(convergencenotattained) {
      #define N (9 + Plasticity_MaxNbOfHardeningVariables + Plasticity_MaxNbOfCriteria)
      double residu[N+3] ;
      double matrix[N*N] ;
      #undef N
      int nmat = 6 + nhardv + ncrit ;
      
      Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm(plasty,stress,hardv,zeta,matrix) ;

      Plasticity_ResidusOfGenericReturnMappingAlgorithm(plasty,stress,hardv,stress_t,hardv_n,zeta,residu+3) ;

      #if DEBUGGENERICRETURNMAPPING
      {
        printf("\n") ;
        printf("=========\n") ;
        printf("iter = %d\n",iter) ;
        printf("=========\n") ;
      }
      #endif

      #if DEBUGGENERICRETURNMAPPING
      {
        int i ;
        
        printf("\n") ;
        printf("jacobian matrix of the residus:\n") ;
        Math_PrintMatrix(matrix,nmat) ;
        
        printf("\n") ;
        printf("residus:\n") ;
        
        for(i = 0 ; i < nmat ; i++) {
          printf("res%d: %e\n",i,residu[3+i]) ;
        }
      }
      #endif
      
      Math_SolveByGaussElimination(matrix,residu+3,nmat) ;
      
      {
        int i ;
        
        for(i = 0 ; i < 6 ; i++) {
          residu[i] = residu[i+3] ;
        }
        
        Elasticity_Convert6TermStressVectorIntoStressTensor(residu) ;
      }
      
      {
        int i ;
        
        for(i = 0 ; i < 9 ; i++) {
          stress[i] += residu[i] ;
        }
    
        for(i = 0 ; i < nhardv ; i++) {
          hardv[i] += residu[9+i] ;
        }
        
        for(i = 0 ; i < ncrit ; i++) {
          zeta[i] += residu[9+nhardv+i] ;
        }
      }
      
      /* Convergence checks */
      {
        double tol = 1.e-6 ;
        int i ;
        
        convergencenotattained = 0 ;
        
        #if 0
        for(i = 0 ; i < 9 ; i++) {
          if(fabs(residu[i]) > tol*fabs(dstress)) {
            convergencenotattained = 1 ;
            break ;
          }
        }
    
        for(i = 0 ; i < nhardv ; i++) {
          if(fabs(residu[9+i]) > tol*fabs(dhardv[i])) {
            convergencenotattained = 1 ;
            break ;
          }
        }
        #endif
        
        for(i = 0 ; i < ncrit ; i++) {
          if(fabs(residu[9+nhardv+i]) > tol*fabs(zeta[i])) {
            convergencenotattained = 1 ;
          }
        }
      }
      
      #if DEBUGGENERICRETURNMAPPING
      {
        int i ;
        
        printf("\n") ;
        printf("solution:\n") ;
        
        printf("stress tensor:\n") ;
        Math_PrintStressTensor(stress) ;
        
        printf("hardening and non-hardening variables:\n") ;
        for(i = 0 ; i < nthardv ; i++) {
          printf("hardv[%d] = %e\n",i,hardv[i]) ;
        }
        
        printf("zeta variables:\n") ;
        for(i = 0 ; i < ncrit ; i++) {
          printf("zeta = %e\n",zeta[i]) ;
        }
      }
      #endif

      if(iter++ > niter) {
        Message_FatalError("Plasticity_GenericReturnMapping: no convergence") ;
      }
    }
  
    /*
      Plastic strains and plastic mulipliers
    */
    {
      Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
      double* flow = flowrules(plasty,stress,hardv) ;
      int k ;
      
      #define FLOW(k)   (flow + nflows*(k))
      for(k = 0 ; k < ncrit ; k++) {
        double* flowk = FLOW(k) ;
        double lambda = Math_Max(zeta[k],0) ;
  
        if(lambda > 0) {
          int    i ;
    
          for(i = 0 ; i < 9 ; i++) {
            strain_p[i] += lambda*flowk[i] ;
          }
        }
    
        Plasticity_GetPlasticMultiplier(plasty)[k] = lambda ;
      }
      
      //Plasticity_FreeBufferFrom(plasty,flow) ;
      #undef FLOW
    }
  
    {
      Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
      double* yield = yieldfunction(plasty,stress_t,hardv_n) ;
      int k ;
    
      for(k = 0 ; k < ncrit ; k++) {
        crit[k] = yield[k] ;
      }
    }
  }
    
  Plasticity_FreeBuffer(plasty) ;
  
  return(crit) ;
}
#endif



#if 1
double* (Plasticity_GenericTangentStiffnessTensor)(Plasticity_t* plasty,const double* stress,const double* hardv,const double* zeta)
/** Compute the consistent tangent stiffness tensor in a generic way
 * 
 *  Inputs are: 
 *    - the current stress tensor
 *    - the current hardening variables
 *    - the current zeta variables (plastic multipliers or yield values)
 * 
 *  On outputs the following values in plasty are modified:
 *    - the consistent tangent stiffness tensor
 * 
 *  Return the value of the yield functions. 
 **/
{
  double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
  double* crit   = Plasticity_GetCriterionValue(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }
  
  if(ncrit > Plasticity_MaxNbOfCriteria) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }

  {
    #define N (6 + Plasticity_MaxNbOfHardeningVariables + Plasticity_MaxNbOfCriteria)
    double matrix[N*N] ;
    #undef N
    int nmat = 6 + nhardv + ncrit ;
      
    Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm(plasty,stress,hardv,zeta,matrix) ;

    #if 0
    {
      printf("\n") ;
      printf("jacobian matrix of the residus:\n") ;
      Math_PrintMatrix(matrix,nmat) ;
    }
    #endif
      
    Math_InvertMatrix(matrix,nmat) ;
      
    #define MATRIX(i,j)  matrix[(i)*nmat+(j)]
    {
      int i ;
        
      for(i = 0 ; i < 6 ; i++) {
        int j ;
          
        for(j = 0 ; j < 6 ; j++) {
          c[6*i + j] = MATRIX(i,j) ;
        }
      }
    }
    #undef MATRIX
      
    Elasticity_Convert6x6MatrixIntoStiffnessMatrix(c) ;
  }
  
  {
    Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
    double* yield = yieldfunction(plasty,stress,hardv) ;
    int k ;
    
    for(k = 0 ; k < ncrit ; k++) {
      crit[k] = yield[k] ;
    }
  }
  
  Plasticity_FreeBuffer(plasty) ;
  
  return(crit) ;
}
#endif



#if 0
double* (Plasticity_GenericTangentStiffnessTensor1)(Plasticity_t* plasty,const double* stress,const double* hardv,const double* plambda)
/** Compute the consistent tangent stiffness tensor in a generic way
 * 
 *  Inputs are: 
 *    - the current stress tensor
 *    - the current hardening variables
 *    - the current plastic multiplier
 *    - a small stress increment used in the numerical derivative
 *    - small increment of hardening variables 
 *    - the number of hardening variables
 * 
 *  On outputs the following values in plasty are modified:
 *    - dfsds = derivative of the yield function wrt stresses
 *    - dgsds = derivative of the potential function wrt stresses
 *    - hm    = hardening modulus
 *    - c     = consistent tangent stiffness tensor
 * 
 *  Return the value of the yield function. 
 **/
{
  Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  
  double* crit   = Plasticity_GetCriterionValue(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  int nflows  = 9 + nhardv ;
  double lambda = (plambda) ? plambda[0] : 0 ;
  
  if(!yieldfunction) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }
  
  if(!flowrules) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }
  
  if(nthardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }


  #define DYIELD_STRESS          (dyield)
  #define DYIELD_HARDV           (dyield+9)
  #define FLOWSTRAIN             (flow)
  #define FLOWHARDV              (flow+9)
  #define DFLOWSTRAIN_STRESS     (dflow)
  #define DFLOWSTRAIN_HARDV      (dflow+81)
  #define DFLOWHARDV_STRESS      (dflow+9*nflows)
  #define DFLOWHARDV_HARDV       (dflow+9*nflows+9)
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    double* dyield = Plasticity_DerivativeOfYieldFunction(yieldfunction,plasty,stress,hardv) ;
    double* flow = flowrules(plasty,stress,hardv) ;
    double* dflow = Plasticity_DerivativeOfFlowRules(flowrules,plasty,stress,hardv) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
    
    {
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        dfsds[i] = DYIELD_STRESS[i] ;
        dgsds[i] = FLOWSTRAIN[i] ;
      }
      
      hm[0] = 0 ;
          
      for(i = 0 ; i < nhardv ; i++) {
        hm[0] -= DYIELD_HARDV[i] * FLOWHARDV[i] ;
      }
    }
    
    #define INVC(i,j)    invc[(i)*9+(j)]
    if(lambda > 0) {
      double* invc = Elasticity_InvertStiffnessMatrix(c) ;

      {
        int    i ;
    
        for(i = 0 ; i < 9 ; i++) {
          int j ;
      
          for(j = 0 ; j < 9 ; j++) {
            INVC(i,j) += lambda * DFLOWSTRAIN_STRESS[i*9+j] ;
          }
        }
      }
        
      #define DFLOWHARDV_STRESS1   (dflowhardvdstress1)
      #define NMAX   Plasticity_MaxNbOfHardeningVariables
      {
        double dflowhardvdstress1[9*NMAX] ;
        double flowhardv1[NMAX] ;
        
        /* Compute flowhardv1 and dflowhardvdstress1 */
        #define D(i,j)        d[(i)*nhardv+(j)]
        {
          double d[NMAX*NMAX] ;
          int k ;
          
          for(k = 0 ; k < nhardv ; k++) {
            int l ;
            
            for(l = 0 ; l < nhardv ; l++) {
              D(k,l) = - lambda * DFLOWHARDV_HARDV[k*nflows+l] ;
            }

            D(k,k) += 1 ;
          }
          
          #define INVD(i,j)        invd[(i)*nhardv+(j)]
          {
            double* invd = Math_InvertMatrix(d,nhardv) ;

            for(k = 0 ; k < nhardv ; k++) {
              int l ;
            
              flowhardv1[k] = 0 ;
            
              for(l = 0 ; l < nhardv ; l++) {
                flowhardv1[k] += INVD(k,l) * FLOWHARDV[l] ;
              }
            }
          
            
            for(k = 0 ; k < nhardv ; k++) {
              int i ;
              
              for(i = 0 ; i < 9 ; i++) {
                int l ;
              
                DFLOWHARDV_STRESS1[k*9+i] = 0 ;
              
                for(l = 0 ; l < nhardv ; l++) {
                  DFLOWHARDV_STRESS1[k*9+i] += INVD(k,l) * DFLOWHARDV_STRESS[l*nflows+i] ;
                }
              }
            }
          }
          #undef INVD
        }
        #undef D

        {
          double lambda2 = lambda * lambda ;
          int k ;
              
          for(k = 0 ; k < nhardv ; k++) {
            int i ;
          
            for(i = 0 ; i < 9 ; i++) {
              int j ;
      
              for(j = 0 ; j < 9 ; j++) {
                INVC(i,j) += lambda2 * DFLOWSTRAIN_HARDV[k*9+i] * DFLOWHARDV_STRESS1[k*9+j] ;
              }
            }
          }
        }

        {
          int k ;
          
          for(k = 0 ; k < nhardv ; k++) {
            int    i ;
    
            for(i = 0 ; i < 9 ; i++) {
              dfsds[i] += lambda * DYIELD_HARDV[k] * DFLOWHARDV_STRESS1[k*9+i] ;
              dgsds[i] += lambda * DFLOWSTRAIN_HARDV[k*9+i] * flowhardv1[k] ;
            }
          }
      
          hm[0] = 0 ;
          
          for(k = 0 ; k < nhardv ; k++) {
            hm[0] -= DYIELD_HARDV[k] * flowhardv1[k] ;
          }
        }
      }
      #undef NMAX
      #undef DFLOWHARDV_STRESS1

      Elasticity_InvertStiffnessMatrix(invc) ;
    }
    #undef INVC
    
    Plasticity_FreeBufferFrom(plasty,dyield) ;
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  #undef FLOWSTRAIN
  #undef FLOWHARDV
  #undef DFLOWSTRAIN_STRESS
  #undef DFLOWSTRAIN_HARDV
  #undef DFLOWHARDV_STRESS
  #undef DFLOWHARDV_HARDV
  #undef DYIELD_HARDV
  #undef DYIELD_STRESS
  
  {
    double* yield  = yieldfunction(plasty,stress,hardv) ;
    
    crit[0] = yield[0] ;
  }
    
  Plasticity_FreeBuffer(plasty) ;
  
  return(crit) ;
}
#endif





/* Herefater is to test the generic approaches */
#if 0
#include <time.h>
#include "Math_.h"

Plasticity_t* Plasticity_DataElasticity(int argc,char** argv)
{
  Plasticity_t* plasty = Plasticity_Create() ;
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double* c = Elasticity_GetStiffnessTensor(elasty) ;
  
  srand(time(NULL)) ;
    
  #if 1
  {
    int rmax = RAND_MAX ;
    double bulk  = (argc > 1) ? (double) atof(argv[1]) : (double) rand() ;
    //double shear = (argc > 2) ? (double) atof(argv[2]) : (double) rand() ;
    //double poisson = (3*bulk - 2*shear)/(6*bulk + 2*shear) ;
    //double young = 9*bulk*shear/(3*bulk + shear) ;
    double poisson = (argc > 2) ? (double) atof(argv[2]) : (double) 0.5*rand()/rmax ;
    double young = 3*bulk*(1 - 2*poisson) ;
    double shear = young/(2 + 2*poisson) ;
    
    printf("Bulk's modulus: %g\n",bulk) ;
    printf("Shear modulus: %g\n",shear) ;
    printf("Young's modulus: %g\n",young) ;
    printf("Poisson's ratio: %g\n",poisson) ;
  
    Elasticity_SetToIsotropy(elasty) ;
    Elasticity_SetParameters(elasty,young,poisson) ;
    
    Elasticity_UpdateElasticTensors(elasty) ;
    
    //Elasticity_ComputeStiffnessTensor(elasty,c) ;
  }
  #endif
  
  #if 0
  {
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
    
    Elasticity_UpdateComplianceTensor(elasty) ;
    
    {
      double bulk = 0 ;
      int i ;

      #define C(i,j)  c[(i)*9 + (j)]
      for(i = 0 ; i < 9 ; i++) {
        int j ;
        
        for(j = 0 ; j < 3 ; j++) {
          bulk += C(3*i+i,3*j+j) ;
        }
      }
      #undef C
      
      Elasticity_GetBulkModulus(elasty) = bulk/9 ;
    }
  }
  #endif
  
  printf("original elastic c:\n") ;
  Math_PrintStiffnessTensor(c) ;
  
  return(plasty) ;
}


double* Plasticity_DataDruckerPrager(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double af = 44*M_PI/180. ;
    double ad = 17*M_PI/180. ;
    double coh = 0 ;
    
    printf("\n") ;
    printf("Drucker-Prager model:\n") ;
    printf("---------------------\n") ;
  
    Plasticity_SetTo(plasty,DruckerPrager) ;
    Plasticity_SetParameters(plasty,af,ad,coh,NULL) ;
  
    {
      int rmax = RAND_MAX ;
      double ff  = 6.*sin(af)/(3. - sin(af)) ;
      double dd  = 6.*sin(ad)/(3. - sin(ad)) ;
      double cc  = 6.*cos(af)/(3. - sin(af))*coh ;
      double k   = 1.e9 ;
      double poisson = 0.2 ;
      double eps = 1.e-3 ;
      double p   = cc/ff + k*eps ;
      double q0  = (ff*p - cc)*1.5*(3 - 6*poisson)/(1 + poisson)/(ff*dd) ;
      double q   = 0.3*q0 ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = ((double) (2*rand()/rmax - 1))*q ;
      double sxz = ((double) (2*rand()/rmax - 1))*q ;
      double syz = ((double) (2*rand()/rmax - 1))*q ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataCamclay(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double lambda = 0.037 ;
    double M = 1.2 ;
    double pc0 = 20.e3 ;
    double phi0 = 0.25 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.004 ;
    
    printf("\n") ;
    printf("CamClay model:\n") ;
    printf("--------------\n") ;
  
    Plasticity_SetTo(plasty,CamClay) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0) ;
  
    {
      int rmax = RAND_MAX / 2 ;
      double p = - 0.9 * pc0 ;
      double q = M*sqrt(-p*(p + pc0)) ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = ((double) (rand() - rmax))/rmax*q ;
      double sxz = ((double) (rand() - rmax))/rmax*q ;
      double syz = ((double) (rand() - rmax))/rmax*q ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataCamclayOffset(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double lambda = 0.037 ;
    double M = 1.2 ;
    double pc0 = 20.e3 ;
    double phi0 = 0.25 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.004 ;
  

    printf("\n") ;
    printf("CamClayOffset model:\n") ;
    printf("--------------------\n") ;
    
    Plasticity_SetTo(plasty,CamClayOffset) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0) ;
  
    {
      int rmax = RAND_MAX / 2 ;
      double ps = 0 ;
      double p = - 0.9 * pc0 ;
      double q = M*sqrt(-(p - ps)*(p + pc0)) ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = ((double) (rand() - rmax))/rmax*q ;
      double sxz = ((double) (rand() - rmax))/rmax*q ;
      double syz = ((double) (rand() - rmax))/rmax*q ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
      hardv[1] = ps ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataBBM(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double lambda = 0.037 ;
    double M = 1.2 ;
    double pc0 = 20.e3 ;
    double phi0 = 0.25 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.004 ;
    double p_ref = 0.01e6 ;
    double coh = 0.8 ;
    Curves_t* curves = Curves_Create(1) ;
    int n = Curves_ReadCurves(curves,"Curves = LC   pc = Range{x1 = 0 , x2 = 1.e6, n = 200} lc = Expressions(1){l0 = 0.037 ; k = 0.004 ; beta = 20.e-6 ; r = 0.75 ; lc = (l0 - k)/(l0*((1-r)*exp(-beta*pc) + r) - k)}") ;
    Curve_t* curvelc   = Curves_FindCurve(curves,"lc") ;
        
    
    printf("\n") ;
    printf("BBM model:\n") ;
    printf("----------\n") ;
    
    Plasticity_SetTo(plasty,BBM) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0,coh,p_ref,curvelc) ;
  
    {
      double lnpc0 = Plasticity_GetHardeningVariable(plasty)[0] ;
      double s     = 100.e3 ;
      double ps    = coh*s ;
      double lc_s  = Curve_ComputeValue(curvelc,s) ;
      double lnp_r = log(p_ref) ;
      double lnpc  = lnp_r + lc_s * (lnpc0 - lnp_r) ;
      double pc    = exp(lnpc) ;
      int rmax = RAND_MAX / 2 ;
      double p = - 0.9 * pc ;
      double q = M*sqrt(-(p - ps)*(p + pc)) ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = ((double) (rand() - rmax))/rmax*q ;
      double sxz = ((double) (rand() - rmax))/rmax*q ;
      double syz = ((double) (rand() - rmax))/rmax*q ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
      hardv[1] = s ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataBExM(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double lambda = 0.016 ;
    double M = 1. ;
    double pc0 = 450.e3 ;
    double eM0 = 0.932 ;
    double em0 = 0.285 ;
    double kappa = 0.045 ;
    double p_ref = 42.e3 ;
    double coh = 0.0073 ;
    double alpha = 2 ;
    //double s_d = 3660.e3 ;
    double s_d = 60.e3 ;
    double s_i = s_d + 25.e3 ;
    double kappa_m = 0.01 ;
    Curves_t* curves = Curves_Create(5) ;
    int n1 = Curves_ReadCurves(curves,"Curves = sl   pc = Range{x1 = 0.01e3 , x2 = 1500.e3, n = 200}  sl = Expressions(1){sl = 0.08+0.92*(1.22-0.22*log(pc/1000))}") ;
    int n2 = Curves_ReadCurves(curves,"Curves = lc   pc = Range{x1 = 0 , x2 = 200.e3, n = 200}  lc = Expressions(1){l0 = 0.14 ; k = 0.015 ; beta = 5.44e-6 ; r = 0.564 ; lc = (l0 - k)/(l0*((1-r)*exp(-beta*pc) + r) - k)}") ;
    int n4 = Curves_ReadCurves(curves,"Curves = fi   p_ratio = Range{x1 = 0 , x2 = 1, n = 100}  fi = Expressions(1){fi = 0.01}") ;
    int n5 = Curves_ReadCurves(curves,"Curves = fd   p_ratio = Range{x1 = 0 , x2 = 1, n = 100}  fd = Expressions(1){fd = 0.41}") ;
    Curve_t* curvelc   = Curves_FindCurve(curves,"lc") ;
    Curve_t* curvefi   = Curves_FindCurve(curves,"fi") ;
    Curve_t* curvefd   = Curves_FindCurve(curves,"fd") ;
    Curve_t* curvesl   = Curves_FindCurve(curves,"sl") ;
        
    
    printf("\n") ;
    printf("BExM model:\n") ;
    printf("-----------\n") ;
    
    Plasticity_SetTo(plasty,BExM) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,eM0,em0,coh,p_ref,alpha,s_i,s_d,kappa_m,curvelc,curvefi,curvefd,curvesl) ;
  
    {
      double lnpc0 = Plasticity_GetHardeningVariable(plasty)[0] ;
      double s     = 100.e3 ;
      double ps    = coh*s ;
      double lc_s  = Curve_ComputeValue(curvelc,s) ;
      double sl_s  = Curve_ComputeValue(curvesl,s) ;
      double lnp_r = log(p_ref) ;
      double lnpc  = lnp_r + lc_s * (lnpc0 - lnp_r) ;
      double pc    = exp(lnpc) ;
      int rmax = RAND_MAX / 2 ;
      double p = - 0.9 * pc ;
      double q = M*sqrt(-(p - ps)*(p + pc)) ;
      //double q = 0 ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = ((double) (rand() - rmax))/rmax*q ;
      double sxz = ((double) (rand() - rmax))/rmax*q ;
      double syz = ((double) (rand() - rmax))/rmax*q ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
      hardv[1] = Plasticity_GetHardeningVariable(plasty)[1] ;
      hardv[2] = Plasticity_GetHardeningVariable(plasty)[2] ;
      hardv[3] = s ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataNSFS(Plasticity_t* plasty,double* stress,double* hardv)
{
  
  {
    double lambda = 0.03 ;
    double M = 1.3 ;
    double pc0 = 15.e6 ;
    double phi0 = 0.15 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.005 ;
    double p_ref = 1.e6 ;
    double coh = 0.8 ;
    double refstrainrate = 1.e-10 ;
    double viscexp = 0 ; // 0.03
    Curves_t* curves = Curves_Create(4) ;
    int n = Curves_ReadCurves(curves,"Curves = wrc   pc = Range{x1 = 0 , x2 = 1.e9, n = 1000}  sl = Expressions(1){p0 = 50e6 ; m = 0.5 ; sl = (1 + (pc/p0)**(1/(1-m)))**(-m)}  kl = Expressions(1){p0 = 50e6 ; m = 0.5 ; A = 2 ; kl = ((1 + (pc/p0)**(1/(1-m)))**(-m))**A} kg = Expressions(1){kg = 1}") ;
    int n1 = Curves_ReadCurves(curves,"Curves = lc   pc = Range{x1 = 0 , x2 = 1.e9, n = 1000}  lc = Expressions(1){l0 = 0.03 ; k = 0.005 ; beta = 0.217e-6 ; r = 0.316 ; lc = (l0 - k)/(l0*((1-r)*exp(-beta*pc) + r) - k)}") ;
    Curve_t* curvesl   = Curves_FindCurve(curves,"sl") ;
    Curve_t* curvelc   = Curves_FindCurve(curves,"lc") ;
    
    
    printf("\n") ;
    printf("NSFS model:\n") ;
    printf("-----------\n") ;
    
    Plasticity_SetTo(plasty,NSFS) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0,coh,p_ref,refstrainrate,viscexp,curvelc,curvesl) ;
  
    {
      double s = 0 ;
      int rmax = RAND_MAX / 2 ;
      double p = - 13.e6 ;
      double q = 16.7e6 ;
      double sxx = p - q/3 ;
      double syy = sxx ;
      double szz = p + 2*q/3 ;
      double sxy = 0 ; //((double) (rand() - rmax))/rmax*q ;
      double sxz = 0 ; //((double) (rand() - rmax))/rmax*q ;
      double syz = 0 ; //((double) (rand() - rmax))/rmax*q ;
      double dt  = 4.e3 ;
      double strainrate = 0 ;
    
      stress[0] = sxx ;
      stress[1] = sxy ;
      stress[2] = sxz ;
      stress[3] = sxy ;
      stress[4] = syy ;
      stress[5] = syz ;
      stress[6] = sxz ;
      stress[7] = syz ;
      stress[8] = szz ;
      
      hardv[0] = Plasticity_GetHardeningVariable(plasty)[0] ;
      hardv[1] = s ;
      hardv[2] = dt ;
      hardv[3] = strainrate ;
    }
  }
  
  return(stress) ;
}



int Plasticity_Tests(Plasticity_t* plasty,const double* stress_t,const double* hardv_n,double* (*ReturnMapping)(Plasticity_t*,double*,double*,double*),double* (*ComputeTangentStiffnessTensor)(Plasticity_t*,const double*,const double*,const double*))
{
  int ncrit  = Plasticity_GetNbOfCriteria(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nnhardv = Plasticity_GetNbOfNonHardeningVariables(plasty) ;
  int nthardv = nhardv + nnhardv ;
  double stress[9] ;
  double hardv[Plasticity_MaxNbOfHardeningVariables] ;
  
  
  printf("Return mapping algorithm:\n") ;
  printf("------------------------\n") ;
  
  {
    
    printf("Trial stress tensor:\n") ;
    Math_PrintStressTensor(stress_t) ;
    
    printf("Hardening and non-hardening variables:\n") ;
    Math_PrintVector(hardv_n,nthardv) ;
    
    printf("yield functions:\n") ;
    {
      double* yield = Plasticity_YieldFunction(plasty,stress_t,hardv_n) ;
      Math_PrintVector(yield,ncrit) ;
    }
    
    {
      double strain_p[9] = {0,0,0,0,0,0,0,0,0} ;
      
      {
        int i ;
        
        for(i = 0 ; i < 9 ; i++) stress[i] = stress_t[i] ;
        
        for(i = 0 ; i < nthardv ; i++) hardv[i] = hardv_n[i] ;
      }
      
      {
        double* yield = ReturnMapping(plasty,stress,strain_p,hardv) ;
        double* lambda = Plasticity_GetPlasticMultiplier(plasty) ;
        
        printf("stress tensor:\n") ;
        Math_PrintStressTensor(stress) ;
    
        printf("Hardening and non-hardening variables:\n") ;
        Math_PrintVector(hardv,nthardv) ;
      
        printf("Plastic strain:\n") ;
        Math_PrintStressTensor(strain_p) ;
      
        printf("yield criteria:\n") ;
        Math_PrintVector(yield,ncrit) ;
      
        printf("Plastic multipliers:\n") ;
        Math_PrintVector(lambda,ncrit) ;
      }
    }
  }
  
  
  printf("Consistent matrix:\n") ;
  printf("-----------------\n") ;
  
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
    double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
    double* hm     = Plasticity_GetHardeningModulus(plasty) ;
    double* dlambda = Plasticity_GetPlasticMultiplier(plasty) ;
    
    {
      double* crit = ComputeTangentStiffnessTensor(plasty,stress,hardv,dlambda) ;
      
      printf("tangent stiffness tensor:\n") ;
      Math_PrintStiffnessTensor(c) ;
      
      printf("yield function gradient:\n") ;
      Math_PrintStressTensor(dfsds) ;
      printf("potential function gradient:\n") ;
      Math_PrintStressTensor(dgsds) ;
      
      printf("hardening moduli:\n") ;
      Math_PrintMatrix(hm,ncrit) ;
    }
  }
  
  return(0) ;
}


/*
 * Compilation: 
 * without petsc and mpi
 * g++ -gdwarf-2 -g3  -L/home/dangla/Documents/Softwares/bil/bil-master/lib -Wl,-rpath=/home/dangla/Documents/Softwares/bil/bil-master/lib -lbil-3.0.0-Debug -o out -lgfortran
 * with petsc and mpi
 * g++ -gdwarf-2 -g3  -L/home/dangla/Documents/Softwares/bil/bil-master/lib -Wl,-rpath=/home/dangla/Documents/Softwares/bil/bil-master/lib -lbil-3.0.0-Debug -L/usr/lib/x86_64-linux-gnu -lpetsc -lmpi -o out -lgfortran
*/
#include "Session.h"

int main(int argc, char** argv)
{
  Session_Open() ;
  
  {
    Plasticity_t* plasty = Plasticity_DataElasticity(argc,argv) ;
    double stress[9] ;
    double hardv[Plasticity_MaxNbOfHardeningVariables] ;
    
    Plasticity_DataDruckerPrager(plasty,stress,hardv) ;
    //Plasticity_DataCamclay(plasty,stress,hardv) ;
    //Plasticity_DataCamclayOffset(plasty,stress,hardv) ;
    //Plasticity_DataBBM(plasty,stress,hardv) ;
    //Plasticity_DataBExM(plasty,stress,hardv) ;
    //Plasticity_DataNSFS(plasty,stress,hardv) ;
    
    {
      
      printf("\n") ;
      printf("native approach\n") ;
      printf("===============\n") ;
      
      Plasticity_Tests(plasty,stress,hardv,Plasticity_GetReturnMapping(plasty),Plasticity_GetComputeTangentStiffnessTensor(plasty)) ;
    
      printf("\n") ;
      printf("generic approach\n") ;
      printf("================\n") ;
      
      Plasticity_Tests(plasty,stress,hardv,Plasticity_GenericReturnMapping,Plasticity_GenericTangentStiffnessTensor) ;
    }
  }
  
  Session_Close() ;
  return(0) ;
}
#endif
