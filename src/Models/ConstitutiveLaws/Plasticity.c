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


static double* (Plasticity_DerivativeOfFlowRules)(Plasticity_t*,Plasticity_FlowRules_t*,const double*,const double*) ;

static double* (Plasticity_DerivativeOfYieldFunction)(Plasticity_t*,Plasticity_YieldFunction_t*,const double*,const double*) ;

static double* (Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm)(Plasticity_t*,const double*,const double*,const double,double*) ;

static double* (Plasticity_ResidusOfGenericReturnMappingAlgorithm)(Plasticity_t*,const double*,const double*,const double*,const double*,const double,double*) ;


/* Drucker-Prager */
#include "PlasticityModels/Plasticity_DruckerPrager.c"


/* Cam-Clay */
#include "PlasticityModels/Plasticity_CamClay.c"


/* Cam-Clay with tensile strength (offset) */
#include "PlasticityModels/Plasticity_CamClayOffset.c"


/* Barcelona Basic model */
#include "PlasticityModels/Plasticity_BBM.c"


/* NSFS model */
#include "PlasticityModels/Plasticity_NSFS.c"


/* Asymmetric Cam-Clay model */
#include "PlasticityModels/Plasticity_ACC.c"
#undef GSLLIB
#ifdef GSLLIB
#include <gsl/gsl_linalg.h>
#include "PlasticityModels/Plasticity_ACCBraun.c"
#endif




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
    double* c = (double*) Mry_New(double[9]) ;
    
    Plasticity_GetYieldFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the potential function gradient */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Plasticity_GetPotentialFunctionGradient(plasty) = c ;
  }
  
  /* Allocation of space for the hardening variable */
  {
    double* c = (double*) Mry_New(double[2*Plasticity_MaxNbOfHardeningVariables]) ;
    
    Plasticity_GetHardeningVariable(plasty) = c ;
    Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) = c + Plasticity_MaxNbOfHardeningVariables ;
  }
  
  /* Allocation of space for the hardening modulus */
  {
    double* c = (double*) Mry_New(double) ;
    
    Plasticity_GetHardeningModulus(plasty) = c ;
  }
  
  /* Allocation of space for Fji*Cijkl */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Plasticity_GetFjiCijkl(plasty) = c ;
  }
  
  /* Allocation of space for Cijkl*Glk */
  {
    double* c = (double*) Mry_New(double[9]) ;
    
    Plasticity_GetCijklGlk(plasty) = c ;
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
  
  return(plasty) ;
}



void  (Plasticity_Delete)(void* self)
{
  Plasticity_t* plasty = (Plasticity_t*) self ;
  
  {
    char* name = Plasticity_GetCodeNameOfModel(plasty) ;
    free(name) ;
  }
  
  {
    double* c = Plasticity_GetYieldFunctionGradient(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetPotentialFunctionGradient(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetHardeningVariable(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetHardeningModulus(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetFjiCijkl(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetCijklGlk(plasty) ;
    free(c) ;
  }
  
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    free(c) ;
  }
  
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    
    Elasticity_Delete(elasty) ;
    free(elasty) ;
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



void Plasticity_Initialize(Plasticity_t* plasty)
{
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTDruckerPrager ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMDruckerPrager ;
    Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFDruckerPrager ;
    Plasticity_GetFlowRules(plasty)                     = Plasticity_FRDruckerPrager ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPDruckerPrager ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 0 ;
    
  } else if(Plasticity_IsCamClay(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTCamClay ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMCamClay ;
    Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFCamClay ;
    Plasticity_GetFlowRules(plasty)                     = Plasticity_FRCamClay ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPCamClay ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 1 ;
    
  } else if(Plasticity_IsCamClayOffset(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTCamClayOffset ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMCamClayOffset ;
    Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFCamClayOffset ;
    Plasticity_GetFlowRules(plasty)                     = Plasticity_FRCamClayOffset ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPCamClayOffset ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 2 ;
    
  } else if(Plasticity_IsBBM(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTBBM ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMBBM ;
    Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFBBM ;
    Plasticity_GetFlowRules(plasty)                     = Plasticity_FRBBM ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPBBM ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 2 ;
    
  } else if(Plasticity_IsACC(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTACC ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMACC ;
    Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFACC ;
    Plasticity_GetFlowRules(plasty)                     = Plasticity_FRACC ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPACC ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 1 ;
    
  #ifdef GSLLIB
  } else if(Plasticity_IsACCBraun(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTACCBraun ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMACCBraun ;
    //Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFACCBraun ;
    //Plasticity_GetFlowRules(plasty)                     = Plasticity_FRACCBraun ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPACCBraun ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 1 ;
  #endif
    
  } else if(Plasticity_IsNSFS(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTNSFSHao ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMNSFSHao ;
    //Plasticity_GetYieldFunction(plasty)                 = Plasticity_YFNSFSHao ;
    //Plasticity_GetFlowRules(plasty)                     = Plasticity_FRNSFSHao ;
    Plasticity_GetSetParameters(plasty)                 = Plasticity_SPNSFSHao ;
    Plasticity_GetNbOfHardeningVariables(plasty)        = 2 ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  
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





void Plasticity_CopyTangentStiffnessTensor(Plasticity_t* plasty,double* c)
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



double Plasticity_UpdateElastoplasticTensor(Plasticity_t* plasty,double* c)
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
        
      Plasticity_GetFjiCijklGlk(plasty) = fcg ;
      
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
        printf("hm + dF(j,i) * C(i,j,k,l) * dG(l,k) = %e\n",det) ;
        printf("\n") ;
        
        Exception_BackupAndTerminate ;
        return(-1) ;
      }
        
      for(i = 0 ; i < 9 ; i++) {
        int j ;
          
        for(j = 0 ; j < 9 ; j++) {
          CEP(i,j) = CEL(i,j) - cg[i]*fc[j]*det ;
        }
      }
  }
  
  return(det) ;
#undef CEP
#undef CEL
}




void Plasticity_PrintTangentStiffnessTensor(Plasticity_t* plasty)
/** Print the 4th rank tangent elastoplastic tensor.
 **/
{
  double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
  
  printf("\n") ;
  printf("4th rank elastoplastic tensor:\n") ;
  
  Math_PrintStiffnessTensor(c) ;
}



double* (Plasticity_DerivativeOfYieldFunction)(Plasticity_t* plasty,double (*yieldfunction)(Plasticity_t*,const double*,const double*),const double* stress,const double* hardv)
/** Return the derivative of the yield function w.r.t. the stresses 
 *  and the hardening variables. 
 *  The output is an array of 9+nhardv components. The 9 first components 
 *  are the derivatives of the yield function w.r.t. the stresses and the
 *  following components are the derivative w.r.t. the hardening variables.
 */
{
  double dstress = Plasticity_GetTypicalSmallIncrementOfStress(plasty) ;
  double* dhardv = Plasticity_GetTypicalSmallIncrementOfHardeningVariable(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nflows = 9 + nhardv ;
  size_t SizeNeeded = nflows*(sizeof(double)) ;
  double* dyield = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  if(dstress == 0) {
    arret("Plasticity_DerivativeOfYieldFunction") ;
  }
  
  {
    double yield0 = yieldfunction(plasty,stress,hardv) ;
    double stress1[9] ;
    double hardv1[Plasticity_MaxNbOfHardeningVariables] ;
    int i ;
    
    for(i = 0 ; i < nhardv ; i++) {
      hardv1[i] = hardv[i] ;
    }
    
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] = stress[i] ;
    }
    
    /* Derivatives with respect to stresses */
    for(i = 0 ; i < 9 ; i++) {
      stress1[i] += dstress ;
      
      {
        double yield1 = yieldfunction(plasty,stress1,hardv) ;
    
        dyield[i] = (yield1 - yield0)/dstress ;
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
        double yield1 = yieldfunction(plasty,stress,hardv1) ;
    
        dyield[9+i] = (yield1 - yield0)/dhardv[i] ;
      }
      
      hardv1[i] = hardv[i] ;
    }
  }
  
  return(dyield) ;
}



double* (Plasticity_DerivativeOfFlowRules)(Plasticity_t* plasty,double* (*flowrules)(Plasticity_t*,const double*,const double*),const double* stress,const double* hardv)
/** Return the derivative of the flow rules w.r.t. the stresses 
 *  and the hardening variables.
 *  The outputs is a (9+nhardv)x(9+nhardv) matrix. The components of 
 *  this matrix are stored in the following order:
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
  int nflows = 9 + nhardv ;
  size_t SizeNeeded = nflows*nflows*(sizeof(double)) ;
  double* dflow = (double*) Plasticity_AllocateInBuffer(plasty,SizeNeeded) ;
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
  
  if(dstress == 0) {
    arret("Plasticity_DerivativeOfFlowRules") ;
  }
  
  #define DFLOWSTRAINDSTRESS(i,j)   (dflow)[(i)*9 + (j)]
  #define DFLOWSTRAINDHARDV(i,j)    (dflow+81)[(j)*9 + (i)]
  #define DFLOWHARDVDSTRESS(i,j)    (dflow+(9+(i))*nflows)[j]
  #define DFLOWHARDVDHARDV(i,j)     (dflow+(9+(i))*nflows+9)[j]
  {
    double* flow0 = flowrules(plasty,stress,hardv) ;
    double stress1[9] ;
    double hardv1[Plasticity_MaxNbOfHardeningVariables] ;
    int i ;
    
    for(i = 0 ; i < nhardv ; i++) {
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
        int j ;
    
        for(j = 0 ; j < 9 ; j++) {
          DFLOWSTRAINDSTRESS(j,i) = (flow1[j] - flow0[j])/dstress ;
        }
    
        for(j = 0 ; j < nhardv ; j++) {
          DFLOWHARDVDSTRESS(j,i) = (flow1[9+j] - flow0[9+j])/dstress ;
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
        int j ;
    
        for(j = 0 ; j < 9 ; j++) {
          DFLOWSTRAINDHARDV(j,i) = (flow1[j] - flow0[j])/dhardv[i] ;
        }
    
        for(j = 0 ; j < nhardv ; j++) {
          DFLOWHARDVDHARDV(j,i) = (flow1[9+j] - flow0[9+j])/dhardv[i] ;
        }
        
        Plasticity_FreeBufferFrom(plasty,flow1) ;
      }
      
      hardv1[i] = hardv[i] ;
    }
    
    Plasticity_FreeBufferFrom(plasty,flow0) ;
  }
  #undef DFLOWSTRAINDSTRESS
  #undef DFLOWSTRAINDHARDV
  #undef DFLOWHARDVDSTRESS
  #undef DFLOWHARDVDHARDV
  
  return(dflow) ;
}



#if 1
double (Plasticity_GenericTangentStiffnessTensor)(Plasticity_t* plasty,const double* stress,const double* hardv,const double dlambda)
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
  
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nflows  = 9 + nhardv ;
  
  if(!yieldfunction) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }
  
  if(!flowrules) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericTangentStiffnessTensor") ;
  }


  #define DYIELDDSTRESS(i)          (dyield)[i]
  #define DYIELDDHARDV(i)           (dyield+9)[i]
  #define FLOWSTRAIN(i)             (flow)[i]
  #define FLOWHARDV(i)              (flow+9)[i]
  #define DFLOWSTRAINDSTRESS(i,j)   (dflow)[(i)*9 + (j)]
  #define DFLOWSTRAINDHARDV(i,j)    (dflow+81+(j)*9)[i]
  #define DFLOWHARDVDSTRESS(i,j)    (dflow+(9+(i))*nflows)[j]
  #define DFLOWHARDVDHARDV(i,j)     (dflow+(9+(i))*nflows+9)[j]
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    double* dyield = Plasticity_DerivativeOfYieldFunction(plasty,yieldfunction,stress,hardv) ;
    double* flow = flowrules(plasty,stress,hardv) ;
    double* dflow = Plasticity_DerivativeOfFlowRules(plasty,flowrules,stress,hardv) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
    
    {
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        dfsds[i] = DYIELDDSTRESS(i) ;
        dgsds[i] = FLOWSTRAIN(i) ;
      }
      
      hm[0] = 0 ;
          
      for(i = 0 ; i < nhardv ; i++) {
        hm[0] -= DYIELDDHARDV(i) * FLOWHARDV(i) ;
      }
    }
    
    #define INVC(i,j)    invc[(i)*9+(j)]
    if(dlambda > 0) {
      double* invc = Elasticity_InvertStiffnessMatrix(c) ;

      {
        int    i ;
    
        for(i = 0 ; i < 9 ; i++) {
          int j ;
      
          for(j = 0 ; j < 9 ; j++) {
            INVC(i,j) += dlambda * DFLOWSTRAINDSTRESS(i,j) ;
          }
        }
      }
        
      #define DFLOWHARDVDSTRESS1(i,j)   (dflowhardvdstress1+(i)*9)[j]
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
              D(k,l) = - dlambda * DFLOWHARDVDHARDV(k,l) ;
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
                flowhardv1[k] += INVD(k,l) * FLOWHARDV(l) ;
              }
            }
          
            
            for(k = 0 ; k < nhardv ; k++) {
              int i ;
              
              for(i = 0 ; i < 9 ; i++) {
                int l ;
              
                DFLOWHARDVDSTRESS1(k,i) = 0 ;
              
                for(l = 0 ; l < nhardv ; l++) {
                  DFLOWHARDVDSTRESS1(k,i) += INVD(k,l) * DFLOWHARDVDSTRESS(l,i) ;
                }
              }
            }
          }
          #undef INVD
        }
        #undef D

        {
          double dlambda2 = dlambda * dlambda ;
          int k ;
              
          for(k = 0 ; k < nhardv ; k++) {
            int i ;
          
            for(i = 0 ; i < 9 ; i++) {
              int j ;
      
              for(j = 0 ; j < 9 ; j++) {
                INVC(i,j) += dlambda2 * DFLOWSTRAINDHARDV(i,k)*DFLOWHARDVDSTRESS1(k,j) ;
              }
            }
          }
        }

        {
          int k ;
          
          for(k = 0 ; k < nhardv ; k++) {
            int    i ;
    
            for(i = 0 ; i < 9 ; i++) {
              dfsds[i] += dlambda * DYIELDDHARDV(k) * DFLOWHARDVDSTRESS1(k,i) ;
              dgsds[i] += dlambda * DFLOWSTRAINDHARDV(i,k) * flowhardv1[k] ;
            }
          }
      
          hm[0] = 0 ;
          
          for(k = 0 ; k < nhardv ; k++) {
            hm[0] -= DYIELDDHARDV(k) * flowhardv1[k] ;
          }
        }
      }
      #undef NMAX
      #undef DFLOWHARDVDSTRESS1

      Elasticity_InvertStiffnessMatrix(invc) ;
    }
    #undef INVC
    
    Plasticity_FreeBufferFrom(plasty,dyield) ;
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  #undef FLOWSTRAIN
  #undef FLOWHARDV
  #undef DFLOWSTRAINDHARDV
  #undef DFLOWHARDVDSTRESS
  #undef DFLOWHARDVDHARDV
  #undef DFLOWSTRAINDSTRESS
  #undef DYIELDDHARDV
  #undef DYIELDDSTRESS
  
  {
    double crit    = yieldfunction(plasty,stress,hardv) ;
    
    return(crit) ;
  }
}
#endif



#if 1
double* (Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm)(Plasticity_t* plasty,const double* stress,const double* hardv,const double dlambda,double* matrix)
/** Compute the jacobian matrix of the function residus 
 *  used for the return mapping algorithm. 
 *  This a (6+nhardv+1)x(6+nhardv+1) matrix associated to
 *    - the plastic strains
 *    - the hardening variables
 *    - the plastic multiplier
 *  Return the pointer to the matrix
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  int nflows  = 9 + nhardv ;
  
  if(!yieldfunction) {
    arret("Plasticity_GenericMatrixOfReturnMappingAlgorithm") ;
  }
  
  if(!flowrules) {
    arret("Plasticity_GenericMatrixOfReturnMappingAlgorithm") ;
  }
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericMatrixOfReturnMappingAlgorithm") ;
  }
  
  
  {
    #define DYIELDDSTRESS             (dyield)
    #define DYIELDDHARDV              (dyield+9)
    #define FLOWSTRAIN                (flow)
    #define FLOWHARDV                 (flow+9)
    #define DFLOWSTRAINDSTRESS        (dflow)
    #define DFLOWSTRAINDHARDV(j)      (dflow+81+(j)*9)
    #define DFLOWHARDVDSTRESS(i)      (dflow+(9+(i))*nflows)
    #define DFLOWHARDVDHARDV(i)       (dflow+(9+(i))*nflows+9)
    #define MATRIX(i,j)               matrix[(i)*nmat+(j)]
    {
      double* dyield = Plasticity_DerivativeOfYieldFunction(plasty,yieldfunction,stress,hardv) ;
      double* flow = flowrules(plasty,stress,hardv) ;
      double* dflow = Plasticity_DerivativeOfFlowRules(plasty,flowrules,stress,hardv) ;
      int nmat = 6 + nhardv + 1 ;
      
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
      
        Elasticity_ConvertStiffnessMatrixInto6x6Matrix(DFLOWSTRAINDSTRESS) ;
        
        for(i = 0 ; i < 6 ; i++) {
          int j ;
          
          for(j = 0 ; j < 6 ; j++) {
            MATRIX(i,j) = invc66[6*i+j] + dlambda * DFLOWSTRAINDSTRESS[6*i+j] ;
          }
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          int i ;
          double* dflowk = DFLOWSTRAINDHARDV(k) ;
          
          Elasticity_ConvertStressTensorInto6TermStressVector(dflowk) ;
          
          for(i = 0 ; i < 6 ; i++) {
            MATRIX(i,6+k) = dlambda * dflowk[i] ;
          }
        }
      }
      
      {
        int i ;
          
        Elasticity_ConvertStressTensorInto6TermStressVector(flow) ;
        
        for(i = 0 ; i < 6 ; i++) {
          MATRIX(i,6+nhardv) = FLOWSTRAIN[i] ;
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          int j ;
          double* dflowk = DFLOWHARDVDSTRESS(k) ;
          
          Elasticity_ConvertStressTensorInto6TermStressVector(dflowk) ;
          
          for(j = 0 ; j < 6 ; j++) {
            MATRIX(6+k,j) = - dlambda * dflowk[j] ;
          }
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          double* dflowk = DFLOWHARDVDHARDV(k) ;
          int l ;
          
          MATRIX(6+k,6+k) = 1 ;
          
          for(l = 0 ; l < nhardv ; l++) {
            MATRIX(6+k,6+l) -= dlambda * dflowk[l] ;
          }
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          MATRIX(6+k,6+nhardv) = - FLOWHARDV[k] ;
        }
      }
      
      {
        int j ;
        
        Elasticity_ConvertStressTensorInto6TermStressVector(DYIELDDSTRESS) ;
        
        for(j = 0 ; j < 6 ; j++) {
          MATRIX(6+nhardv,j) = DYIELDDSTRESS[j] ;
        }
      }
      
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          MATRIX(6+nhardv,6+k) = DYIELDDHARDV[k] ;
        }
      }
      
      Plasticity_FreeBufferFrom(plasty,dyield) ;
    }
    #undef FLOWSTRAIN
    #undef FLOWHARDV
    #undef DFLOWSTRAINDHARDV
    #undef DFLOWHARDVDSTRESS
    #undef DFLOWHARDVDHARDV
    #undef DFLOWSTRAINDSTRESS
    #undef DYIELDDHARDV
    #undef DYIELDDSTRESS
    #undef MATRIX
  }
  
  return(matrix) ;
}
#endif



#if 1
double* (Plasticity_ResidusOfGenericReturnMappingAlgorithm)(Plasticity_t* plasty,const double* stress,const double* hardv,const double* stress_t,const double* hardv_n,const double dlambda,double* residu)
/** Compute the residu vector used for the return mapping algorithm. 
 *  
 *  On input residu should point to an array of 9+nhardv+1 doubles.
 * 
 *  On output the (6+nhardv+1) first components of the residu vector 
 *  are associated to
 *    - the plastic strain flow rule
 *    - the hardening flow rule
 *    - the yield function
 * 
 *  Return the pointer to the residu
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  
  if(!yieldfunction) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  if(!flowrules) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_ResidusOfGenericReturnMappingAlgorithm") ;
  }
  
  
  {
    double* invc = Elasticity_GetComplianceTensor(elasty) ;
    
    #define FLOWSTRAIN                (flow)
    #define FLOWHARDV                 (flow+9)
    #define RESIDUSTRAIN              (residu)
    #define RESIDUHARDV               (residu+6)
    #define RESIDUYIELD               (residu+6+nhardv)
    {
      double yield = yieldfunction(plasty,stress,hardv) ;
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
          
          RESIDUSTRAIN[i] = dlambda * FLOWSTRAIN[i] ;
          
          for(j = 0 ; j < 9 ; j++) {
            RESIDUSTRAIN[i] += invc[9*i + j] * dstress[j] ;
          }
        }
        
        Elasticity_ConvertStressTensorInto6TermStressVector(RESIDUSTRAIN) ;
      }
      
      /* Equation of the hardening flow rules */
      {
        int k ;
        
        for(k = 0 ; k < nhardv ; k++) {
          RESIDUHARDV[k] = hardv[k] - hardv_n[k] - dlambda * FLOWHARDV[k] ;
        }
      }
      
      /* Equation of the yield criterion */
      {
        RESIDUYIELD[0] = yield ;
      }
      
      /* Change the sign */
      {
        int nmat = 6 + nhardv + 1 ;
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
double (Plasticity_GenericReturnMapping)(Plasticity_t* plasty,double* stress,double* strain_p,double* hardv)
/** Return mapping algorithm in a generic way
 * 
 *  On inputs: 
 *    - "stress" points to the current elastic trial stress tensor
 *    - "hardv" points to the hardening variables obtained at the previous time step
 *    - "strain_p" points to the plastic strains obtained at the previous time step
 * 
 *  On outputs, the following values are modified:
 *    - the stresses ("stress"), 
 *    - the plastic strains ("strain_p"), 
 *    - the hardening variables ("hardv").
 * 
 *  Return the value of the yield function. 
 **/
{
  Plasticity_YieldFunction_t* yieldfunction = Plasticity_GetYieldFunction(plasty) ;
  int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
  double yield = 0 ;
  double dlambda = 0 ;
  
  if(!yieldfunction) {
    arret("Plasticity_GenericReturnMapping") ;
  }
  
  if(nhardv > Plasticity_MaxNbOfHardeningVariables) {
    arret("Plasticity_GenericReturnMapping") ;
  }
  
  yield = yieldfunction(plasty,stress,hardv) ;
  
  if(yield > 0.) {
    double hardv_n[Plasticity_MaxNbOfHardeningVariables] ;
    double stress_t[9] ;
    int iter = 0 ;
    int niter = 20 ;
    int convergencenotattained = 1 ;
  
    {
      int i ;
    
      for(i = 0 ; i < 9 ; i++) {
        stress_t[i] = stress[i] ;
      }
    
      for(i = 0 ; i < nhardv ; i++) {
        hardv_n[i] = hardv[i] ;
      }
    }
    
    while(convergencenotattained) {
      #define N (9 + Plasticity_MaxNbOfHardeningVariables + 1)
      double residu[N+3] ;
      double matrix[N*N] ;
      #undef N
      int nmat = 6 + nhardv + 1 ;
      
      Plasticity_JacobianMatrixOfGenericReturnMappingAlgorithm(plasty,stress,hardv,dlambda,matrix) ;

      Plasticity_ResidusOfGenericReturnMappingAlgorithm(plasty,stress,hardv,stress_t,hardv_n,dlambda,residu+3) ;

      #if 0
      {
        int i ;
        
        printf("\n") ;
        printf("=========\n") ;
        printf("iter = %d\n",iter) ;
        printf("=========\n") ;
        
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
        
        dlambda += residu[9+nhardv] ;
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
        
        {
          if(fabs(residu[9+nhardv]) > tol*fabs(dlambda)) {
            convergencenotattained = 1 ;
          }
        }
      }
      
      #if 0
      {
        int i ;
        
        printf("stress tensor:\n") ;
        Math_PrintStressTensor(stress) ;
        
        for(i = 0 ; i < nhardv ; i++) {
          printf("hardv[%d] = %e\n",i,hardv[i]) ;
        }
        
        printf("lambda = %e\n",dlambda) ;
      }
      #endif

      if(iter++ > niter) {
        Message_FatalError("Plasticity_GenericReturnMapping: no convergence") ;
      }
    }
  }
  
  /*
    Plastic strains
  */
  
  if(dlambda > 0) {
    Plasticity_FlowRules_t* flowrules = Plasticity_GetFlowRules(plasty) ;
  
    if(flowrules) {
      double* flow = flowrules(plasty,stress,hardv) ;
      int    i ;
    
      for(i = 0 ; i < 9 ; i++) {
        strain_p[i] += dlambda*flow[i] ;
      }
      
      Plasticity_FreeBufferFrom(plasty,flow) ;
    } else {
      arret("Plasticity_GenericReturnMapping") ;
    }
  }
  
  /* Plastic muliplier */
  Plasticity_GetPlasticMultiplier(plasty) = dlambda ;
  
  return(yield) ;
}
#endif






#if 0
static int Plasticity_TestMatrix(Plasticity_t*,const double*) ;
static int Plasticity_TestReturnMapping(Plasticity_t*,const double*) ;

#include <time.h>
#include "Math_.h"

Plasticity_t* Plasticity_DataElasticity(int argc,char** argv)
{
  Plasticity_t* plasty = Plasticity_Create() ;
  
  srand(time(NULL)) ;
  srand(rand()) ;
    
  #if 1
  {
    double bulk  = (argc > 1) ? (double) atof(argv[1]) : (double) rand() ;
    double shear = (argc > 2) ? (double) atof(argv[2]) : (double) rand() ;
    double poisson = (3*bulk - 2*shear)/(6*bulk + 2*shear) ;
    double young = 9*bulk*shear/(3*bulk+shear) ;
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
    printf("Bulk's modulus: %g\n",bulk) ;
    printf("Shear modulus: %g\n",shear) ;
    printf("Young's modulus: %g\n",young) ;
    printf("Poisson's ratio: %g\n",poisson) ;
  
    Elasticity_SetToIsotropy(elasty) ;
    Elasticity_SetParameters(elasty,young,poisson) ;
    
    Elasticity_UpdateElasticTensors(elasty) ;
    
    //Elasticity_ComputeStiffnessTensor(elasty,c) ;
  
    printf("original elastic c:\n") ;
    Math_PrintStiffnessTensor(c) ;
  }
  #endif
  
  #if 0
  {
    Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
    double* c = Elasticity_GetStiffnessTensor(elasty) ;
    
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
  
    printf("original elastic c:\n") ;
    Math_PrintStiffnessTensor(c) ;
  }
  #endif
  
  return(plasty) ;
}


double* Plasticity_DataCamclay(Plasticity_t* plasty)
{
  double* stress = (double*) Mry_New(double[9]) ;
  
  {
    double lambda = 0.037 ;
    double M = 1.2 ;
    double pc0 = 20.e3 ;
    double phi0 = 0.25 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.004 ;
  
    Plasticity_SetToCamClay(plasty) ;
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
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataCamclayOffset(Plasticity_t* plasty)
{
  double* stress = (double*) Mry_New(double[9]) ;
  
  {
    double lambda = 0.037 ;
    double M = 1.2 ;
    double pc0 = 20.e3 ;
    double phi0 = 0.25 ;
    double e0 = phi0/(1 - phi0) ;
    double kappa = 0.004 ;
  
    Plasticity_SetToCamClay(plasty) ;
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
      
      Plasticity_GetHardeningVariable(plasty)[1] = ps ;
    }
  }
  
  return(stress) ;
}


double* Plasticity_DataBBM(Plasticity_t* plasty)
{
  double* stress = (double*) Mry_New(double[9]) ;
  
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
    Curve_t* lc   = Curves_FindCurve(curves,"lc") ;
        
    Plasticity_SetToBBM(plasty) ;
    Plasticity_SetParameters(plasty,kappa,lambda,M,pc0,e0,coh,p_ref,lc) ;
  
    {
      double* hardv = Plasticity_GetHardeningVariable(plasty) ;
      double m     = Plasticity_GetSlopeCriticalStateLine(plasty) ;
      double k     = Plasticity_GetSuctionCohesionCoefficient(plasty) ;
      double p_r   = Plasticity_GetReferenceConsolidationPressure(plasty) ;
      Curve_t* lc  = Plasticity_GetLoadingCollapseFactorCurve(plasty) ;
      double lnpc0 = hardv[0] ;
      double s     = 0 ;
      double ps    = k*s ;
      double lc_s  = Curve_ComputeValue(lc,s) ;
      double lnp_r = log(p_r) ;
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
      
      Plasticity_GetHardeningVariable(plasty)[1] = s ;
    }
  }
  
  return(stress) ;
}


int Plasticity_TestMatrix(Plasticity_t* plasty,const double* stress)
{
  
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
    double dlambda = 1 ;
    double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
    double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
    double* hm     = Plasticity_GetHardeningModulus(plasty) ;
    double* hardv =  Plasticity_GetHardeningVariable(plasty) ;
    
    {
      double crit = Plasticity_ComputeTangentStiffnessTensor(plasty,stress,hardv,dlambda) ;
      
      printf("\n") ;
      printf("native approach:\n") ;
      printf("---------------\n") ;
      printf("tangent stiffness tensor:\n") ;
      Math_PrintStiffnessTensor(c) ;
      
      printf("yield function gradient:\n") ;
      Math_PrintStressTensor(dfsds) ;
      printf("potential function gradient:\n") ;
      Math_PrintStressTensor(dgsds) ;
      printf("hm[0] = %e\n",hm[0]) ;
    }
    
    {
      double crit = Plasticity_GenericTangentStiffnessTensor(plasty,stress,hardv,dlambda) ;
    
      printf("\n") ;
      printf("generic approach:\n") ;
      printf("----------------\n") ;
      printf("tangent stiffness tensor:\n") ;
      Math_PrintStiffnessTensor(c) ;
      
      printf("yield function gradient:\n") ;
      Math_PrintStressTensor(dfsds) ;
      printf("potential function gradient:\n") ;
      Math_PrintStressTensor(dgsds) ;
      printf("hm[0] = %e\n",hm[0]) ;
    }
  }
  
  return(0) ;
}



int Plasticity_TestReturnMapping(Plasticity_t* plasty,const double* stress_t)
{
  
  {
    double* hardv_n = Plasticity_GetHardeningVariable(plasty) ;
    
    printf("Trial stress tensor:\n") ;
    Math_PrintStressTensor(stress_t) ;
    printf("hardv = %e\n",hardv_n[0]) ;
    printf("yield function:\n") ;
    printf("Yield criterion = %e\n",Plasticity_YieldFunction(plasty,stress_t,hardv_n)) ;
    
    {
      double stress[9] ;
      double hardv[Plasticity_MaxNbOfHardeningVariables] ;
      double strain_p[9] = {0,0,0,0,0,0,0,0,0} ;
      double yield ;
      
      {
        int i ;
        int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
        
        for(i = 0 ; i < 9 ; i++) stress[i] = stress_t[i] ;
        
        for(i = 0 ; i < nhardv ; i++) hardv[i] = hardv_n[i] ;
      }
      
      printf("\n") ;
      printf("native approach:\n") ;
      printf("---------------\n") ;
      
      yield = Plasticity_ReturnMapping(plasty,stress,strain_p,hardv) ;
      
      printf("stress tensor:\n") ;
      Math_PrintStressTensor(stress) ;
      printf("hardv = %e\n",hardv[0]) ;
      printf("Plastic strain:\n") ;
      Math_PrintStressTensor(strain_p) ;
      printf("Yield criterion = %e\n",yield) ;
      printf("Plastic multiplier = %e\n",Plasticity_GetPlasticMultiplier(plasty)) ;
    }
    
    {
      double stress[9] ;
      double hardv[Plasticity_MaxNbOfHardeningVariables] ;
      double strain_p[9] = {0,0,0,0,0,0,0,0,0} ;
      double yield ;
      
      {
        int i ;
        int nhardv = Plasticity_GetNbOfHardeningVariables(plasty) ;
        
        for(i = 0 ; i < 9 ; i++) stress[i] = stress_t[i] ;
        
        for(i = 0 ; i < nhardv ; i++) hardv[i] = hardv_n[i] ;
      }
    
      printf("\n") ;
      printf("generic approach:\n") ;
      printf("----------------\n") ;
      
      yield = Plasticity_GenericReturnMapping(plasty,stress,strain_p,hardv) ;
      
      printf("stress tensor:\n") ;
      Math_PrintStressTensor(stress) ;
      printf("hardv = %e\n",hardv[0]) ;
      printf("Plastic strain:\n") ;
      Math_PrintStressTensor(strain_p) ;
      printf("Yield criterion = %e\n",yield) ;
      printf("Plastic multiplier = %e\n",Plasticity_GetPlasticMultiplier(plasty)) ;
    }
  }
  
  return(0) ;
}


/*
 * Compilation: 
 * g++ -gdwarf-2 -g3  -L/home/dangla/Documents/Softwares/bil/bil-master/lib -Wl,-rpath=/home/dangla/Documents/Softwares/bil/bil-master/lib -lbil-2.8.8-Debug -o out -lgfortran
*/
#include "Session.h"

int main(int argc, char** argv)
{
  Session_Open() ;
  
  {
    Plasticity_t* plasty = Plasticity_DataElasticity(argc,argv) ;
    //double* stress = Plasticity_DataCamclay(plasty) ;
    //double* stress = Plasticity_DataCamclayOffset(plasty) ;
    double* stress = Plasticity_DataBBM(plasty) ;
  
    printf("Test on the consistent matrix\n") ;
    printf("=============================\n") ;
    Plasticity_TestMatrix(plasty,stress) ;
  
    printf("Test on the return mapping algorithm\n") ;
    printf("====================================\n") ;
    Plasticity_TestReturnMapping(plasty,stress) ;
  }
  
  Session_Close() ;
  return(0) ;
}
#endif
