#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "Exception.h"
#include "Mry.h"
#include "Message.h"
#include "DataFile.h"
#include "Tools/Math.h"
#include "Plasticity.h"


/* Drucker-Prager */
static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_ComputeFunctionGradientsDruckerPrager ;
static Plasticity_ReturnMapping_t                    Plasticity_ReturnMappingDruckerPrager ;


/* Cam-Clay */
static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_ComputeFunctionGradientsCamClay ;
static Plasticity_ReturnMapping_t                    Plasticity_ReturnMappingCamClay ;
static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_ComputeFunctionGradientsCamClayEp ;
static Plasticity_ReturnMapping_t                    Plasticity_ReturnMappingCamClayEp ;
static Plasticity_ComputeTangentStiffnessTensor_t    Plasticity_ComputeFunctionGradientsCamClayOffset ;
static Plasticity_ReturnMapping_t                    Plasticity_ReturnMappingCamClayOffset ;




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
    double* c = (double*) Mry_New(double[Plasticity_MaxNbOfHardeningVariables]) ;
    
    Plasticity_GetHardeningVariable(plasty) = c ;
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
  
  /* Allocation of space for the parameters */
  {
    double* c = (double*) Mry_New(double[Plasticity_MaxNbOfParameters]) ;
    
    Plasticity_GetParameter(plasty) = c ;
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
    double* c = Plasticity_GetParameter(plasty) ;
    free(c) ;
  }
}



void Plasticity_Initialize(Plasticity_t* plasty)
{
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_ComputeFunctionGradientsDruckerPrager ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_ReturnMappingDruckerPrager ;
    
  } else if(Plasticity_IsCamClay(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_ComputeFunctionGradientsCamClay ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_ReturnMappingCamClay ;
    
  } else if(Plasticity_IsCamClayEp(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_ComputeFunctionGradientsCamClayEp ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_ReturnMappingCamClayEp ;
    
  } else if(Plasticity_IsCamClayOffset(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_ComputeFunctionGradientsCamClayOffset ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_ReturnMappingCamClayOffset ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  
}



void Plasticity_SetParameter(Plasticity_t* plasty,const char* str,double v)
{
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    if(0) {
    } else if(!strcmp(str,"friction angle")) {
      Plasticity_GetFrictionAngle(plasty)  = v ;
    } else if(!strcmp(str,"dilatancy angle")) {
      Plasticity_GetDilatancyAngle(plasty) = v ;
    } else if(!strcmp(str,"cohesion")) {
      Plasticity_GetCohesion(plasty)       = v ;
      Plasticity_GetHardeningVariable(plasty)[0] = 0 ;
    }

  } else if(Plasticity_IsCamClay(plasty)) {
    if(0) {
    } else if(!strcmp(str,"slope swelling line")) {
      Plasticity_GetSlopeSwellingLine(plasty)               = v ;
    } else if(!strcmp(str,"slope virgin consolidation line")) {
      Plasticity_GetSlopeVirginConsolidationLine(plasty)    = v ;
    } else if(!strcmp(str,"slope critical state line")) {
      Plasticity_GetSlopeCriticalStateLine(plasty)          = v ;
    } else if(!strcmp(str,"initial preconsolidation pressure")) {
      Plasticity_GetInitialPreconsolidationPressure(plasty) = v ;
      Plasticity_GetHardeningVariable(plasty)[0] = v ;
    }

  } else if(Plasticity_IsCamClayEp(plasty)) {
    if(0) {
    } else if(!strcmp(str,"slope swelling line")) {
      Plasticity_GetSlopeSwellingLine(plasty)               = v ;
    } else if(!strcmp(str,"slope virgin consolidation line")) {
      Plasticity_GetSlopeVirginConsolidationLine(plasty)    = v ;
    } else if(!strcmp(str,"slope critical state line")) {
      Plasticity_GetSlopeCriticalStateLine(plasty)          = v ;
    } else if(!strcmp(str,"initial preconsolidation pressure")) {
      Plasticity_GetInitialPreconsolidationPressure(plasty) = v ;
      Plasticity_GetHardeningVariable(plasty)[0] = 0 ;
    }

  } else if(Plasticity_IsCamClayOffset(plasty)) {
    if(0) {
    } else if(!strcmp(str,"slope swelling line")) {
      Plasticity_GetSlopeSwellingLine(plasty)               = v ;
    } else if(!strcmp(str,"slope virgin consolidation line")) {
      Plasticity_GetSlopeVirginConsolidationLine(plasty)    = v ;
    } else if(!strcmp(str,"slope critical state line")) {
      Plasticity_GetSlopeCriticalStateLine(plasty)          = v ;
    } else if(!strcmp(str,"initial preconsolidation pressure")) {
      Plasticity_GetInitialPreconsolidationPressure(plasty) = v ;
      Plasticity_GetHardeningVariable(plasty)[0] = v ;
    }
    
  } else {
    Message_RuntimeError("Not known") ;
  }

}



void Plasticity_SetParameters(Plasticity_t* plasty,...)
{
  va_list args ;
  
  va_start(args,plasty) ;
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    Plasticity_GetFrictionAngle(plasty)            = va_arg(args,double) ;
    Plasticity_GetDilatancyAngle(plasty)           = va_arg(args,double) ;
    Plasticity_GetCohesion(plasty)                 = va_arg(args,double) ;
    
    Plasticity_GetHardeningVariable(plasty)[0] = 0 ;
    
  } else if(Plasticity_IsCamClay(plasty)) {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialVoidRatio(plasty)                = va_arg(args,double) ;
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      
      Plasticity_GetHardeningVariable(plasty)[0] = pc ;
    }
    
  } else if(Plasticity_IsCamClayEp(plasty)) {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialVoidRatio(plasty)                = va_arg(args,double) ;
    
    Plasticity_GetHardeningVariable(plasty)[0] = 0 ;
    
  } else if(Plasticity_IsCamClayOffset(plasty)) {
    Plasticity_GetSlopeSwellingLine(plasty)               = va_arg(args,double) ;
    Plasticity_GetSlopeVirginConsolidationLine(plasty)    = va_arg(args,double) ;
    Plasticity_GetSlopeCriticalStateLine(plasty)          = va_arg(args,double) ;
    Plasticity_GetInitialPreconsolidationPressure(plasty) = va_arg(args,double) ;
    Plasticity_GetInitialVoidRatio(plasty)                = va_arg(args,double) ;
    
    {
      double pc = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
      
      Plasticity_GetHardeningVariable(plasty)[0] = pc ;
    }
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  

  va_end(args) ;
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
      Curves_FreeBuffer(curves) ;
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






#if 0
void Plasticity_CopyElasticTensor(Plasticity_t* plasty,double* c)
/** Copy the 4th rank elastic tensor in c. */
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double* cel = Elasticity_GetStiffnessTensor(elasty) ;
  
  /* Elastic stiffness tensor */
  {
      int i ;
        
      for(i = 0 ; i < 81 ; i++) {
        c[i] = cel[i] ;
      }
  }
}
#endif








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





/* Local functions */
double Plasticity_ComputeFunctionGradientsDruckerPrager(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
/** Drucker-Prager criterion. 
 * 
 *  Inputs are: 
 *  the stresses (sig),
 *  the cumulative plastic shear strain (gam_p = hardv[0], not used here). 
 * 
 *  Parameters are:
 *  the friction angle (af),
 *  the dilatancy angle (ad) and the cohesion.
 * 
 *  On outputs the following values are modified:
 *  dfsds = derivative of the yield function wrt stresses
 *  dgsds = derivative of the potential function wrt stresses
 *  hm    = hardening modulus
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  //double gam_p   = hardv[0] ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,dev[9],devn[9] ;
  double crit ;
  double cc ;
  
  
  /*
    Cohesion
  */
  {
    //double c1 = (gam_p < gam_R) ? 1 - (1 - alpha)*gam_p/gam_R : alpha ;
    double c1 = 1 ;
    
    cc = cc0*c1*c1 ;
  }
  
  /*
    Yield function
  */ 
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q + ff*p - cc ;
  
  /*
    Deviatoric stresses
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dev[i] = sig[i] - p*id[i] ;
    }
    
    if(q > 0.) {
      for(i = 0 ; i < 9 ; i++) {
        devn[i] = 1.5 * dev[i] / q ;
      }
    } else {
      for(i = 0 ; i < 9 ; i++) {
        devn[i] = 0 ;
      }
    }
  }
  
  /*
    Yield function gradient
  */
  if(q > 0.) {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = 1.5*dev[i]/q + id[i]*ff/3. ;
    }
    
  } else {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      dfsds[i] = id[i]*ff/3. ;
    }
  }
  
  /*
    Potential function gradient
  */
  
  /* Elastic case */
  if(crit <= 0.) {
    if(q > 0.) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    } else {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = id[i]*dd/3. ;
      }
    }
  }
  
  /* Plastic case */
  if(crit > 0.) {
    double k   = young/(3 - 6*poisson) ;
    double dmu = young/(1 + poisson) ;
    double mu  = 0.5*dmu ;
    
    /* Smooth flow regime */
    if(q > crit*3*mu/(3*mu+k*ff*dd)) {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = 1.5*dev[i]/q + id[i]*dd/3. ;
      }
      
    /* Flow regime at the notch apex */
    } else {
      double dl = (ff*p - cc)/(k*ff*dd) ;
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        dgsds[i] = dev[i]/(dmu*dl) + id[i]*dd/3. ;
      }
    }
  }
  
  /* Hardening modulus */
  hm[0] = 0. ;
  /*
  if(gam_p < gam_R) {
    hm[0] = -2.*(1.-alpha)/gam_R*(1.-(1.-alpha)*gam_p/gam_R)*cc0 ;
    hm[0] *= sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(dgsds)) ;
  }
  */
  
  Plasticity_GetCriterionValue(plasty) = crit ;
  
  /*
   * Consistent tangent matrix
   * Before updating we modify the elastic matrix according to
   * C0_ijkl += - Dlambda*4*G^2/(2q) * ( 3 * I_ijkl - 1_ij * 1_kl - 2 * devn_ij * devn_kl )
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
     
    if(dlambda > 0 && q > 0.) {
      double bulk    = Elasticity_GetBulkModulus(elasty) ;
      double shear   = Elasticity_GetShearModulus(elasty) ;
      double g0 = shear ;
      double k0 = bulk  ;
      double g1 = g0 * q / (3*g0*dlambda + q) ;
      double k1 = k0 ;
   
      {
        double mu      = g1 ;
        double lame    = k1 - 2./3.*mu ;
        int    i ;

        for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
        for(i = 0 ; i < 3 ; i++) {
          int j ;
      
          #define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
          for(j = 0 ; j < 3 ; j++) {
            C(i,i,j,j) += lame ;
            C(i,j,i,j) += mu ;
            C(i,j,j,i) += mu ;
          }
          #undef C
        }
      }
       
      {
        double a = - 4*g1*g1*dlambda/(3*g1*dlambda - q) ;
        int    i ;
        
        for(i = 0 ; i < 9 ; i++) {
          int j ;
         
          #define C1(i,j)    (c[(i)*9+(j)])
          for(j = 0 ; j < 9 ; j++) {
            C1(i,j) += a*devn[i]*devn[j] ;
          }
          #undef C1
        }
      }
    }
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
   
  return(crit) ;
}



double Plasticity_ReturnMappingDruckerPrager(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Drucker-Prager return mapping.
 * 
 *  Parameters are:
 *  the Young modulus (young),
 *  the Poisson's ratio (poisson),
 *  the friction angle (af), 
 *  the dilatancy angle (ad),
 *  the cohesion (cohesion).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the cumulative plastic shear strain (gam_p = hardv[0]).
 *  the plastic multiplier in plasty.
 * 
 *  Return the value of the yield function. */
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = dmu/2. ;
  double k       = young/(3 - 6*poisson) ;
  double af      = Plasticity_GetFrictionAngle(plasty) ;
  double ad      = Plasticity_GetDilatancyAngle(plasty) ;
  double cohesion = Plasticity_GetCohesion(plasty) ;
  double ff      = 6.*sin(af)/(3. - sin(af)) ;
  double dd      = 6.*sin(ad)/(3. - sin(ad)) ;
  double cc0     = 6.*cos(af)/(3. - sin(af))*cohesion ;
  double gam_p   = hardv[0] ;

  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p_t,q_t,sdev_t[9] ;
  double crit ;
  double dl = 0 ;
  
  /*
    Trial stresses
  */ 
  p_t  = (sig[0] + sig[4] + sig[8])/3. ;
  q_t  = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      sdev_t[i] = sig[i] - p_t*id[i] ;
    }
  }
  
  /*
    Criterion
  */ 
  {
    //double c1 = ((gam_p) < gam_R) ? 1 - (1 - alpha)*(gam_p)/gam_R : alpha ;
    double c1 = 1 ;
    double cc = cc0*c1*c1 ;
    
    crit = q_t + ff*p_t - cc ;
  }
  
  /*
    Return mapping: update plastic strains and stresses
  */
  if(crit > 0.) {
    double deps_p[9] = {0,0,0,0,0,0,0,0,0} ;
    double p = p_t ;
    double q = q_t ;
    double gam_pn = gam_p ;
    double gam_p1 = gam_p ;
    
    /* Smooth flow regime: assuming that q > 0 */
    if(q > 0) {
      double fcrit = crit ;
      int    nf    = 0 ;
      double tol   = 1.e-10 ;
      
      //while(fabs(fcrit) > tol*cc0) {
      while(fabs(fcrit/q) > tol) {
        double dqsdl ;
        double dpsdl ;
        double dccsdl ;
        double c1 ;
        double cc ;
        int    i ;
        
        /* Plastic strain increments */
        for(i = 0 ; i < 9 ; i++) {
          deps_p[i] = dl*(1.5*sdev_t[i]/q_t + id[i]*dd/3.) ;
        }
        
        /* Cumulative plastic shear strain */
        gam_p1 = gam_pn + sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(deps_p)) ;
        
        /* p, q, cc */
        //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
        c1 = 1 ;
        cc = cc0*c1*c1 ;
        q  = q_t - dl*3*mu ;
        p  = p_t - dl*k*dd ;
        
        /* dqsdl, dpsdl, dccsdl */
        dqsdl = -3*mu ;
        dpsdl = -k*dd ;
        dccsdl = 0. ;
        /*
        if(gam_p1 < gam_R) {
          dccsdl = -2*(1 - alpha)/gam_R*c1*cc0 ;
          dccsdl *= 1.5*sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(sdev_t))/q_t ;
        }
        */
        
        /* Criterion */
        fcrit = q + ff*p - cc ;
        
        /* dl */
        {
          double df = dqsdl + ff*dpsdl - dccsdl ;
          
          dl   -= fcrit/df ;
        }
        
        if(nf++ > 20) {
          Message_FatalError("Plasticity_ReturnMappingDruckerPrager: no convergence") ;
        }
      }
      
      /* Stresses */
      {
        int    i ;
        
        for(i = 0 ; i < 9 ; i++) {
          sig[i] = sdev_t[i]*q/q_t + p*id[i] ;
        }
      }
    }
    
    /* Flow regime at the notch apex */
    if(q <= 0.) {
      double c1 ;
      double cc ;
      int    i ;
      
      /* Deviatoric plastic strain increments */
      for(i = 0 ; i < 9 ; i++) deps_p[i]  = sdev_t[i]/dmu ;
        
      /* Cumulative plastic shear strain */
      gam_p1 = gam_pn + sqrt(2*Math_ComputeSecondDeviatoricStressInvariant(deps_p)) ;
      
      /* p, q, cc */
      //c1 = (gam_p1 < gam_R) ? 1 - (1 - alpha)*gam_p1/gam_R : alpha ;
      c1 = 1 ;
      cc = cc0*c1*c1 ;
      p  = cc/ff ;
      q  = 0. ;
      
      /* dl */
      dl   = (ff*p_t - cc)/(k*ff*dd) ;
      
      /* Plastic strain increments and stresses */
      for(i = 0 ; i < 9 ; i++) {
        deps_p[i] = sdev_t[i]/dmu + dl*id[i]*dd/3. ;
        sig[i]    = p*id[i] ;
      }
    }
    
      
    /* Total plastic strains */
    {
      int    i ;
      
      for(i = 0 ; i < 9 ; i++) {
        eps_p[i] += deps_p[i] ;
      }
    }
    hardv[0] = gam_p1 ;
  }
  
  /* Plastic multiplier */
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}





double Plasticity_ComputeFunctionGradientsCamClayEp(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
/** Modified Cam-Clay criterion */
{
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  double e_p     = hardv[0] ;
  double pc      = pc0 * exp(-v*e_p) ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Gradients
    ---------
    dp/dsig_ij = 1/3 delta_ij
    dq/dsig_ij = 3/2 dev_ij/q 
    df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q 
    df/dp      = 2*p + pc
    df/dq      = 2*q/m2
    
    df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij 
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
    
      dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
      dgsds[i] = dfsds[i] ;
    }
  }
  
  /* The hardening modulus */
  {
    //double v = 1./(lambda - kappa) ;
    
    hm[0] = (1 + e0)*v*p*(2*p + pc)*pc ;
  }
  
  /*
   * Tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  
  return(crit) ;
}



double Plasticity_ReturnMappingCamClayEp(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Modified Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Dangla.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (m),
 *  the initial pre-consolidation pressure (pc),
 *  the initial void ratio (e0).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the plastic void ratio (e_p = hardv[0]).
 * 
 *  Return the value of the yield function. 
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = 0.5*dmu ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double phi0    = e0/(1 + e0) ;
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  double e_p     = hardv[0] ;
  double pc      = pc0 * exp(-v*e_p) ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q ;
  double p_t,q_t ;
  double crit ;
  double dl ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + p*(p + pc) = 0
     for p. The other variables (pc,q,dl) are expressed with p.
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.e-8 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      /*
       * Flow rule
       * ---------
       * df/dp      = 2*p + pc
       * df/dq      = (2/m2) q
       * dp/dsig_ij = 1/3 delta_ij
       * dq/dsig_ij = 3/2 dev_ij/q
       * df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij 
       * deps_p     = dl (2*p + pc)
       * deij_p     = dl (3/m2) dev_ij
       */
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double ddlsdp = ((1 - phi0)*kappa/p - dl*(2+dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*ddlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      /* Plastic multiplier (dl):
       * ------------------------
       * deps_e = - (1-phi0) kappa ln(p/p_n) ; 
       * deps   = - (1-phi0) kappa ln(p_t/p_n) ; 
       * deps_p = - (1-phi0) kappa ln(p_t/p) = dl (2*p + pc)
       * Hence 
       * dl = (1-phi0) kappa ln(p/p_t) / (2*p + pc)
       * Pre-consolidation pressure (pc):
       * --------------------------------
       * deps_p = - (1-phi0) kappa ln(p_t/p) = - (1-phi0) (lambda - kappa) ln(pc/pc_n)
       * Hence 
       * ln(pc/pc_n) = - v kappa ln(p/p_t)
       * Deviatoric behavior (q):
       * ------------------------
       * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
       * Hence 
       * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
       * q      = q_t / (1 + 6 mu / m2 dl)
       */
       
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
        Message_FatalError("Plasticity_ReturnMappingCamClayEp: no convergence") ;
      }
    }
  }
  
  /*
    Stresses and plastic strains
  */
  
  {
    double a = 1./(1 + 6*mu/m2*dl) ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev      = a*(sig[i] - p_t*id[i]) ;
      double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
      sig[i]    = p*id[i] + dev ;
      eps_p[i] += dl*dfsds ;
    }
  }
  
  /* Hardening variable */
  {
    double de_p = (1 + e0) * dl * (2*p + pc) ;
    
    hardv[0] += de_p ;
  }
  
  /* Plastic muliplier */
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}





double Plasticity_ComputeFunctionGradientsCamClay(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
/** Modified Cam-Clay criterion */
{
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double pc      = hardv[0] ;
  double m2      = m*m ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
    Gradients
    ---------
    dp/dsig_ij = 1/3 delta_ij
    dq/dsig_ij = 3/2 dev_ij/q 
    df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q 
    df/dp      = 2*p + pc
    df/dq      = 2*q/m2
    
    df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij 
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
    
      dfsds[i] = (2*p + pc)*id[i]/3. + 3./m2*dev ;
      dgsds[i] = dfsds[i] ;
    }
  }
  
  /* The hardening modulus */
  /* H is defined by: df = (df/dsig_ij) dsig_ij - dl H 
   * But df = (df/dsig_ij) dsig_ij + (df/dpc) dpc
   * Hence: H = - (df/dpc) dpc / dl
   * On the other hand 
   * dpc/pc = - v de_p = - (1 + e0) v deps_p = - (1 + e0) v dl (dg/dp)
   * i.e. dpc = dl hpc with hpc = - (1 + e0)*v*pc*(2*p + pc)
   * Hence: H = (1 + e0) v (df/dpc) (dg/dp) pc
   */
  {
    double v = 1./(lambda - kappa) ;
    
    hm[0] = (1 + e0)*v*p*(2*p + pc)*pc ;
  }
  
  /*
   * Tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
     
    if(dlambda > 0) {
      double v        = 1./(lambda - kappa) ;
      double v1       = (1 + e0)*v ;
      double hpc      = - v1*(2*p + pc)*pc ;
      double dhpcdpc  = - v1*2*(p + pc) ;
      double dhpcdp   = - v1*2*pc ;
      double dfdpc    = p ;
      double ddgdpdpc = 1 ;
      double dlambda1 = dlambda / (1 - dlambda*dhpcdpc) ;
      Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
      double bulk    = Elasticity_GetBulkModulus(elasty) ;
      double shear   = Elasticity_GetShearModulus(elasty) ;
      double g0 = shear ;
      double k0 = bulk  ;
      double g1 = g0 / (6*g0*dlambda/m2 + 1) ;
      double k1 = k0 / (2*k0*dlambda*(1 - dlambda1*v1*pc) + 1) ;
   
      {
        double mu      = g1 ;
        double lame    = k1 - 2./3.*mu ;
        int    i ;

        for(i = 0 ; i < 81 ; i++) c[i] = 0. ;
    
        for(i = 0 ; i < 3 ; i++) {
          int j ;
      
          #define C(i,j,k,l)  (c[(((i)*3+(j))*3+(k))*3+(l)])
          for(j = 0 ; j < 3 ; j++) {
            C(i,i,j,j) += lame ;
            C(i,j,i,j) += mu ;
            C(i,j,j,i) += mu ;
          }
          #undef C
        }
      }
      
      {
        int    i ;
    
        for(i = 0 ; i < 9 ; i++) {
          dfsds[i] += dlambda1*dfdpc*dhpcdp*id[i]/3 ;
          dgsds[i] += dlambda1*hpc*ddgdpdpc*id[i]/3 ;
        }
      }
      
      hm[0] /= (1 - dlambda*dhpcdpc) ;
    }
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  
  return(crit) ;
}



double Plasticity_ReturnMappingCamClay(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Modified Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Dangla.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc=hardv[0]),
 *  the initial void ratio (e0).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the pre-consolidation pressure (pc=hardv[0]).
 * 
 *  Return the value of the yield function. 
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = 0.5*dmu ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double pc      = hardv[0] ;
  double phi0    = e0/(1 + e0) ;
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  double p_t,q_t ;
  double dl ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + p*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + p*(p + pc) = 0
     for p. The other variables (pc,q,dl) are expressed with p.
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.e-8 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      /*
       * Flow rule
       * ---------
       * df/dp      = 2*p + pc
       * df/dq      = (2/m2) q
       * dp/dsig_ij = 1/3 delta_ij
       * dq/dsig_ij = 3/2 dev_ij/q
       * df/dsig_ij = 1/3 (2*p + pc) delta_ij + (3/m2) dev_ij 
       * deps_p     = dl (2*p + pc)
       * deij_p     = dl (3/m2) dev_ij
       */
      double dfsdp  = 2*p + pc ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p ;
      double dpcsdp = -v*kappa*pc/p ;
      double ddlsdp = ((1 - phi0)*kappa/p - dl*(2 + dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*ddlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      /* Plastic multiplier (dl):
       * ------------------------
       * deps_e = - (1-phi0) kappa ln(p/p_n) ; 
       * deps   = - (1-phi0) kappa ln(p_t/p_n) ; 
       * deps_p = - (1-phi0) kappa ln(p_t/p) = dl (2*p + pc)
       * Hence 
       * dl = (1-phi0) kappa ln(p/p_t) / (2*p + pc)
       * Pre-consolidation pressure (pc):
       * --------------------------------
       * deps_p = - (1-phi0) kappa ln(p_t/p) = - (1-phi0) (lambda - kappa) ln(pc/pc_n)
       * Hence 
       * ln(pc/pc_n) = - v kappa ln(p/p_t)
       * Deviatoric behavior (q):
       * ------------------------
       * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
       * Hence 
       * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
       * q      = q_t / (1 + 6 mu / m2 dl)
       */
       
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + p*(p + pc) ;
      
      if(nf++ > 20) {
        Message_FatalError("Plasticity_ReturnMappingCamClay: no convergence") ;
      }
    }
  }
  
  /*
    Stresses and plastic strains
  */
  
  {
    double a = 1./(1 + 6*mu/m2*dl) ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev      = a*(sig[i] - p_t*id[i]) ;
      double dfsds    = (2*p + pc)*id[i]/3. + 3./m2*dev ;
    
      sig[i]    = p*id[i] + dev ;
      eps_p[i] += dl*dfsds ;
    }
  }
  
  /* Consolidation pressure */
  hardv[0] = pc ;
  
  /* Plastic muliplier */
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}





double Plasticity_ComputeFunctionGradientsCamClayOffset(Plasticity_t* plasty,const double* sig,const double* hardv,const double dlambda)
/** Modified Cam-Clay with offset criterion */
{
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double* dfsds  = Plasticity_GetYieldFunctionGradient(plasty) ;
  double* dgsds  = Plasticity_GetPotentialFunctionGradient(plasty) ;
  double* hm     = Plasticity_GetHardeningModulus(plasty) ;
  double pc      = hardv[0] ;
  double ps      = hardv[1] ;
  double m2      = m*m ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  
  /* 
     The yield criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + (p - ps)*(p + pc) ;
  
  /*
    Gradients
    ---------
    dp/dsig_ij = 1/3 delta_ij
    dq/dsig_ij = 3/2 dev_ij/q 
    df/dsig_ij = 1/3 (df/dp) delta_ij + 3/2 (df/dq) dev_ij/q 
    df/dp      = 2*p + pc - ps
    df/dq      = 2*q/m2
    
    df/dsig_ij = 1/3 (2*p + pc - ps) delta_ij + (3/m2) dev_ij 
  */
  {
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev = sig[i] - p*id[i] ;
    
      dfsds[i] = (2*p + pc - ps)*id[i]/3. + 3./m2*dev ;
      dgsds[i] = dfsds[i] ;
    }
  }
  
  /* The hardening modulus */
  /* H is defined by: df = (df/dsig_ij) dsig_ij - dl H 
   * But df = (df/dsig_ij) dsig_ij + (df/dpc) dpc
   * Hence: H = - (df/dpc) dpc / dl
   * On the other hand 
   * dpc/pc = - v de_p = - (1 + e0) v deps_p = - (1 + e0) v dl (dg/dp)
   * Hence: H = (1 + e0) v (df/dpc) (dg/dp) pc
   */
  {
    double v = 1./(lambda - kappa) ;
    
    hm[0] = (1 + e0)*v*(p - ps)*(2*p + pc - ps)*pc ;
  }
  
  /*
   * Tangent matrix
   */
  {
    double* c = Plasticity_GetTangentStiffnessTensor(plasty) ;
     
    Plasticity_CopyElasticTensor(plasty,c) ;
       
    Plasticity_UpdateElastoplasticTensor(plasty,c) ;
  }
  
  return(crit) ;
}



double Plasticity_ReturnMappingCamClayOffset(Plasticity_t* plasty,double* sig,double* eps_p,double* hardv)
/** Modified Cam-Clay return mapping.
 *  Algorithm from Borja & Lee 1990 modified by Dangla.
 * 
 *  Inputs are: 
 *  the slope of the swelling line (kappa),
 *  the slope of the virgin consolidation line (lambda),
 *  the shear modulus (mu),
 *  the slope of the critical state line (M),
 *  the pre-consolidation pressure (pc),
 *  the initial void ratio (e0).
 *  The tensile strength offsetting (ps=hardv[1]).
 * 
 *  On outputs, the following values are modified:
 *  the stresses (sig), 
 *  the plastic strains (eps_p), 
 *  the pre-consolidation pressure (pc=hardv[0]).
 * 
 *  Return the value of the yield function. 
 **/
{
  Elasticity_t* elasty = Plasticity_GetElasticity(plasty) ;
  double young   = Elasticity_GetYoungModulus(elasty) ;
  double poisson = Elasticity_GetPoissonRatio(elasty) ;
  double dmu     = young/(1 + poisson) ;
  double mu      = 0.5*dmu ;
  double m       = Plasticity_GetSlopeCriticalStateLine(plasty) ;
  double kappa   = Plasticity_GetSlopeSwellingLine(plasty) ;
  double lambda  = Plasticity_GetSlopeVirginConsolidationLine(plasty) ;
  double e0      = Plasticity_GetInitialVoidRatio(plasty) ;
  //double pc0     = Plasticity_GetInitialPreconsolidationPressure(plasty) ;
  double pc      = hardv[0] ;
  double ps      = hardv[1] ;
  double phi0    = e0/(1 + e0) ;
  double m2      = m*m ;
  double v       = 1./(lambda - kappa) ;
  
  double id[9] = {1,0,0,0,1,0,0,0,1} ;
  double p,q,crit ;
  double p_t,q_t ;
  double dl ;
  
  /* 
     The criterion
  */
  p    = (sig[0] + sig[4] + sig[8])/3. ;
  q    = sqrt(3*Math_ComputeSecondDeviatoricStressInvariant(sig)) ;
  crit = q*q/m2 + (p - ps)*(p + pc) ;
  
  /*
     Closest point projection algorithm.
   * Only one iterative loop is used to solve
                    q*q/m2 + (p - ps)*(p + pc) = 0
     for p. The other variables (pc,q,dl) are expressed with p.
   */
  dl    = 0. ;
  p_t   = p ;
  q_t   = q ;
  
  if(crit > 0.) {
    double pc_n  = pc ;
    double fcrit = crit ;
    int nf    = 0 ;
    double tol = 1.e-8 ;
    
    while(fabs(fcrit) > tol*pc_n*pc_n) {
      /*
       * Flow rule
       * ---------
       * df/dp      = 2*p + pc - ps
       * df/dq      = (2/m2) q
       * dp/dsig_ij = 1/3 delta_ij
       * dq/dsig_ij = 3/2 dev_ij/q
       * df/dsig_ij = 1/3 (2*p + pc - ps) delta_ij + (3/m2) dev_ij 
       * deps_p     = dl  (2*p + pc - ps)
       * deij_p     = dl  (3/m2) dev_ij
       */
      double dfsdp  = 2*p + pc - ps ;
      double dfsdq  = 2*q/m2 ;
      double dfsdpc = p - ps ;
      double dpcsdp = -v*kappa*pc/p ;
      double ddlsdp = ((1 - phi0)*kappa/p - dl*(2 + dpcsdp))/dfsdp ;
      double dqsdp  = -q*6*mu/(m2 + 6*mu*dl)*ddlsdp ;
      double df     = dfsdp + dfsdq*dqsdp + dfsdpc*dpcsdp ;
      
      p     -= fcrit/df ;
      
      /* Variables (pc,dl,q) are explicit in p */
      /* Plastic multiplier (dl):
       * ------------------------
       * deps_e = - (1-phi0) kappa ln(p/p_n) ; 
       * deps   = - (1-phi0) kappa ln(p_t/p_n) ; 
       * deps_p = - (1-phi0) kappa ln(p_t/p) = dl (2*p + pc - ps)
       * Hence 
       * dl = (1-phi0) kappa ln(p/p_t) / (2*p + pc - ps)
       * Pre-consolidation pressure (pc):
       * --------------------------------
       * deps_p = - (1-phi0) kappa ln(p_t/p) = - (1-phi0) (lambda - kappa) ln(pc/pc_n)
       * Hence 
       * ln(pc/pc_n) = - v kappa ln(p/p_t)
       * Deviatoric behavior (q):
       * ------------------------
       * dev_ij = dev_ij_t - 2 mu deij_p = dev_ij_t - 6 mu / m2 dl dev_ij
       * Hence 
       * dev_ij = dev_ij_t / (1 + 6 mu / m2 dl)
       * q      = q_t / (1 + 6 mu / m2 dl)
       */
       
      pc     = pc_n*pow(p/p_t,-v*kappa) ;
      dl     = (1 - phi0)*kappa*log(p/p_t)/(2*p + pc - ps) ;
      q      = q_t*m2/(m2 + 6*mu*dl) ;
      fcrit  = q*q/m2 + (p - ps)*(p + pc) ;
      
      if(nf++ > 20) {
        Message_FatalError("Plasticity_ReturnMappingCamClayOffset: no convergence") ;
      }
    }
  }
  
  /*
    Stresses and plastic strains
  */
  
  {
    double a = 1./(1 + 6*mu/m2*dl) ;
    int    i ;
    
    for(i = 0 ; i < 9 ; i++) {
      double dev      = a*(sig[i] - p_t*id[i]) ;
      double dfsds    = (2*p + pc - ps)*id[i]/3. + 3./m2*dev ;
    
      sig[i]    = p*id[i] + dev ;
      eps_p[i] += dl*dfsds ;
    }
  }
  
  /* Consolidation pressure */
  hardv[0] = pc ;
  
  /* Plastic muliplier */
  Plasticity_GetPlasticMultiplier(plasty) = dl ;
  
  return(crit) ;
}








/* Not used */
#if 0
double Plasticity_ComputeFunctionGradients(Plasticity_t* plasty,const double* stress,const double* hardv,const double dlambda)
{
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    return(Plasticity_ComputeFunctionGradientsDruckerPrager(plasty,stress,hardv,dlambda)) ;
  } else if(Plasticity_IsCamClay(plasty)) {
    return(Plasticity_ComputeFunctionGradientsCamClay(plasty,stress,hardv,dlambda)) ;
  } else {
    Message_RuntimeError("Not known") ;
  }
  
  return(0) ;
}



double Plasticity_ReturnMapping(Plasticity_t* plasty,double* stress,double* strain_p,double* hardv)
{
  
  if(Plasticity_IsDruckerPrager(plasty)) {
    return(Plasticity_ReturnMappingDruckerPrager(plasty,stress,strain_p,hardv)) ;
  } else if(Plasticity_IsCamClay(plasty)) {
    return(Plasticity_ReturnMappingCamClay(plasty,stress,strain_p,hardv)) ;
  } else {
    Message_RuntimeError("Not known") ;
  }
  
  return(0) ;
}
#endif
