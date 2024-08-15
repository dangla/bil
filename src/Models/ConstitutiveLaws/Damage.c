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
#include "Damage.h"


extern Damage_SetModelProp_t Damage_ListOfSetModelProp ;



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
    
    Damage_GetTangentStiffnessTensor(damage) = c ;
  }
  
  /* Allocation of space for the damaged stiffness tensor */
  {
    double* c = (double*) Mry_New(double[81]) ;
    
    Damage_GetDamagedStiffnessTensor(damage) = c ;
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
    double* c = Damage_GetTangentStiffnessTensor(damage) ;
    free(c) ;
  }
  
  {
    double* c = Damage_GetDamagedStiffnessTensor(damage) ;
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
  
  if(Damage_Is(damage,"Mazars")) {
    Damage_GetSetModelProp(damage) = DamageMazars_SetModelProp ;
      
  } else if(Damage_Is(damage,"MarigoJirasek")) {
    Damage_GetSetModelProp(damage) = DamageMarigoJirasek_SetModelProp ;
    
  } else {
    Message_RuntimeError("Not known") ;
  }
  
}





void (DamageModel_Initialize)(Damage_t* damage,const char* codename)
{
  int n_models = Damage_NbOfModels ;
  const char* modelnames[] = {Damage_ListOfNames} ;
  Damage_SetModelProp_t* xDamage_SetModelProp[] = {Damage_ListOfSetModelProp} ;
  int i = 0 ;
  
  while(i < n_models && strcmp(modelnames[i],codename)) i++ ;
    
  if(i < n_models) {
    Damage_CopyCodeName(damage,modelnames[i]) ;
    Damage_GetSetModelProp(damage) = xDamage_SetModelProp[i] ;
    Damage_SetModelProp(damage) ;
    
    return ;
  } else {
    
    Message_Warning("No model named %s",codename) ;
  }

  return ;
}



void Damage_CopyTangentStiffnessTensor(Damage_t* damage,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* ct = Damage_GetTangentStiffnessTensor(damage) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = ct[i] ;
    }
  }
}



void Damage_CopyDamagedStiffnessTensor(Damage_t* damage,double* c)
/** Copy the 4th rank stiffness tensor in c. */
{
  double* cd = Damage_GetDamagedStiffnessTensor(damage) ;

  {
    int i ;
        
    for(i = 0 ; i < 81 ; i++) {
      c[i] = cd[i] ;
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
            
        Message_Direct("\n") ;
        Message_Direct("Damage_UpdateTangentStiffnessTensor:") ;
        Message_Direct("\n") ;
        Message_Direct("Something wrong happened\n") ;
            
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
        
        Exception_BackupAndTerminate ;
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




void Damage_PrintTangentStiffnessTensor(Damage_t* damage)
/** Print the 4th rank tangent damage tensor.
 **/
{
  double* c = Damage_GetTangentStiffnessTensor(damage) ;
  
  printf("\n") ;
  printf("4th rank tangent damage tensor:\n") ;
  
  Math_PrintStiffnessTensor(c) ;
}
