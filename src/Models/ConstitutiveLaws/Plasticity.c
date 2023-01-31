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


/* Drucker-Prager */
#include "PlasticityModels/Plasticity_DruckerPrager.c"


/* Cam-Clay */
#include "PlasticityModels/Plasticity_CamClay.c"


/* Cam-Clay with tensile strength (offset) */
#include "PlasticityModels/Plasticity_CamClayOffset.c"




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
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTDruckerPrager ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMDruckerPrager ;
    
  } else if(Plasticity_IsCamClay(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTCamClay ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMCamClay ;
    
  } else if(Plasticity_IsCamClayOffset(plasty)) {
    Plasticity_GetComputeTangentStiffnessTensor(plasty) = Plasticity_CTCamClayOffset ;
    Plasticity_GetReturnMapping(plasty)                 = Plasticity_RMCamClayOffset ;
    
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

  } else if(Plasticity_IsCamClay(plasty) ||
            Plasticity_IsCamClayOffset(plasty)) {
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
    
  } else if(Plasticity_IsCamClay(plasty) ||
            Plasticity_IsCamClayOffset(plasty)) {
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
