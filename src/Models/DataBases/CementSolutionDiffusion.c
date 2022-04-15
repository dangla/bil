#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Message.h"
#include "Exception.h"
#include "Math_.h"
#include "Temperature.h"
#include "CementSolutionChemistry.h"
#include "CementSolutionDiffusion.h"





static void  (CementSolutionDiffusion_AllocateMemory)     (CementSolutionDiffusion_t*) ;
static void  (CementSolutionDiffusion_Initialize)         (CementSolutionDiffusion_t*) ;
static void  (CementSolutionDiffusion_UpdateElementFluxes)(CementSolutionDiffusion_t*) ;




CementSolutionDiffusion_t* CementSolutionDiffusion_Create(void)
{
  CementSolutionDiffusion_t* csd = (CementSolutionDiffusion_t*) malloc(sizeof(CementSolutionDiffusion_t)) ;
  
  if(!csd) arret("CementSolutionDiffusion_Create") ;
  
  {
    /* Memory allocation */
    CementSolutionDiffusion_AllocateMemory(csd) ;
  
    /* Initialize concentrations to zero and other related variables */
    CementSolutionDiffusion_Initialize(csd) ;
  }
  
  return(csd) ;
}



void CementSolutionDiffusion_AllocateMemory(CementSolutionDiffusion_t* csd)
{
  
  
  /* Allocation of space for the temperature */
  {
    Temperature_t* temp = Temperature_Create() ;
    
    CementSolutionDiffusion_GetTemperature(csd) = temp ;
  }
  
  /* Allocation of space for the diffusion coefficients */
  {
    size_t sz = CementSolutionDiffusion_NbOfSpecies*sizeof(double) ;
    double* d = (double*) malloc(sz) ;
    
    if(!d) arret("CementSolutionDiffusion_AllocateMemory(10)") ;
    
    CementSolutionDiffusion_GetDiffusionCoefficient(csd) = d ;
  }
  
  
  /* Allocation of space for the gradients */
  {
    size_t sz = CementSolutionDiffusion_NbOfSpecies*sizeof(double) ;
    double* grd = (double*) malloc(sz) ;
    
    if(!grd) arret("CementSolutionDiffusion_AllocateMemory(11)") ;
    
    CementSolutionDiffusion_GetGradient(csd) = grd ;
  }
  
  
  /* Allocation of space for the potentials */
  {
    int n = CementSolutionDiffusion_MaxNbOfPotentialVectors ;
    size_t sz = n*CementSolutionDiffusion_NbOfSpecies*sizeof(double) ;
    double* pot = (double*) malloc(sz) ;
    
    if(!pot) arret("CementSolutionDiffusion_AllocateMemory(11)") ;
    
    CementSolutionDiffusion_GetPotential(csd) = pot ;
  }
  
  
  /* Allocation of space for the pointer to potentials */
  {
    int n = CementSolutionDiffusion_MaxNbOfPotentialVectors ;
    size_t sz = n*sizeof(double*) ;
    double** ppot = (double**) malloc(sz) ;
    
    if(!ppot) arret("CementSolutionDiffusion_AllocateMemory(11)") ;
    
    CementSolutionDiffusion_GetPointerToPotentials(csd) = ppot ;
    
    {
      double* pot = CementSolutionDiffusion_GetPotential(csd) ;
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        ppot[i] = pot + i*CementSolutionDiffusion_NbOfSpecies ;
      }
    }
  }
  
  
  /* Allocation of space for the concentration fluxes */
  {
    size_t sz = CementSolutionDiffusion_NbOfSpecies*sizeof(double) ;
    double* flx = (double*) malloc(sz) ;
    
    if(!flx) arret("CementSolutionDiffusion_AllocateMemory(11)") ;
    
    CementSolutionDiffusion_GetFlux(csd) = flx ;
  }
  
  
  /* Allocation of space for the element fluxes */
  {
    size_t sz = CementSolutionDiffusion_NbOfElementFluxes*sizeof(double) ;
    double* efx = (double*) malloc(sz) ;
    
    if(!efx) arret("CementSolutionDiffusion_AllocateMemory(6)") ;
    
    CementSolutionDiffusion_GetElementFlux(csd) = efx ;
  }
}




/* Intern functions */

#include "DiffusionCoefficientOfMoleculeInWater.h"

/* Shorthands of macros */
#define Diff(CPD)             CementSolutionDiffusion_GetDiffusionCoefficientOf(csd,CPD)


void CementSolutionDiffusion_Initialize(CementSolutionDiffusion_t* csd)
{
  /* Diffusion coefficients */
  {
    int     n = CementSolutionDiffusion_NbOfSpecies ;
    double* v = CementSolutionDiffusion_GetDiffusionCoefficient(csd) ;
    double TK = CementSolutionDiffusion_GetRoomTemperature(csd) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  
#define D(A)         DiffusionCoefficientOfMoleculeInWater(A,TK)
    Diff(OH)       = D(OH) ;
    Diff(H)        = D(H) ;
  
    Diff(Ca)       = D(Ca) ;
    Diff(CaOH)     = D(CaOH) ;
    Diff(CaO2H2)   = D(CaO2H2) ;
  
    Diff(H4SiO4)   = D(H4SiO4) ;
    Diff(H3SiO4)   = D(H3SiO4) ;
    Diff(H2SiO4)   = D(H2SiO4) ;
  
    Diff(Na)       = D(Na) ;
    Diff(NaOH)     = D(NaOH) ;
  
    Diff(K)        = D(K) ;
    Diff(KOH)      = D(KOH) ;
  
    Diff(Cl)       = D(Cl) ;
  
    Diff(H2CO3)    = D(H2CO3) ;
    Diff(HCO3)     = D(HCO3) ;
    Diff(CO3)      = D(CO3) ;
    Diff(CO2)      = D(CO2) ;
  
    Diff(H2SO4)    = D(H2SO4) ;
    Diff(HSO4)     = D(HSO4) ;
    Diff(SO4)      = D(SO4) ;
  
    Diff(Al)       = D(Al) ;
    Diff(AlO4H4)   = D(AlO4H4) ;
  
    Diff(NaHCO3)   = D(NaHCO3) ;
    Diff(NaCO3)    = D(NaCO3) ;
  
    Diff(CaHCO3)   = D(CaHCO3) ;
    Diff(CaCO3)    = D(CaCO3) ;
  
    Diff(CaH2SiO4) = D(CaH2SiO4) ;
    Diff(CaH3SiO4) = D(CaH3SiO4) ;

    Diff(CaHSO4)   = D(CaHSO4) ;
    Diff(CaSO4)    = D(CaSO4) ;

#undef D
  }
  
  /* Gradients */
  {
    int     n = CementSolutionDiffusion_NbOfSpecies ;
    double* v = CementSolutionDiffusion_GetGradient(csd) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
  
  /* Fluxes */
  {
    int     n = CementSolutionDiffusion_NbOfSpecies ;
    double* v = CementSolutionDiffusion_GetFlux(csd) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
  
  /* Element fluxes */
  {
    int     n = CementSolutionDiffusion_NbOfElementFluxes ;
    double* v = CementSolutionDiffusion_GetElementFlux(csd) ;
  
    {
      int i ;
      
      for(i = 0 ; i < n ; i++) {
        v[i] = 0. ;
      }
    }
  }
}




void CementSolutionDiffusion_ComputeFluxes(CementSolutionDiffusion_t* csd)
{
  /* The fluxes */
  {
    double* grad = CementSolutionDiffusion_GetGradient(csd) ;
    double* flux = CementSolutionDiffusion_GetFlux(csd) ;
    double* diff = CementSolutionDiffusion_GetDiffusionCoefficient(csd) ;
    double* z    = CementSolutionChemistry_GetValence() ;
    int n = CementSolutionDiffusion_NbOfSpecies ;
    double  current = 0. ;
    int i ;

    for(i = 0 ; i < n ; i++) {
      flux[i]  = - diff[i] * grad[i] ;
      current += z[i]*flux[i] ;
    }
    
    CementSolutionDiffusion_GetIonCurrent(csd) = current ;
  }
  
  CementSolutionDiffusion_UpdateElementFluxes(csd) ;
}




/* Shorthands of macros */
#define Flux(CPD)             CementSolutionDiffusion_GetFluxOf(csd,CPD)

void CementSolutionDiffusion_UpdateElementFluxes(CementSolutionDiffusion_t* csd)
/** Update the element concentrations **/
{
  double w_ca       = Flux(Ca) ;
  double w_caoh     = Flux(CaOH) ;
  double w_caoh2    = Flux(CaO2H2) ;
  
  double w_h4sio4   = Flux(H4SiO4) ;
  double w_h3sio4   = Flux(H3SiO4) ;
  double w_h2sio4   = Flux(H2SiO4) ;
  
  double w_na       = Flux(Na) ;
  double w_naoh     = Flux(NaOH) ;
  
  double w_k        = Flux(K) ;
  double w_koh      = Flux(KOH) ;
  
  double w_co2      = Flux(CO2) ;
  double w_h2co3    = Flux(H2CO3) ;
  double w_hco3     = Flux(HCO3) ;
  double w_co3      = Flux(CO3) ;
  
  double w_h2so4    = Flux(H2SO4) ;
  double w_hso4     = Flux(HSO4) ;
  double w_so4      = Flux(SO4) ;
  
  double w_al       = Flux(Al) ;
  double w_alo4h4   = Flux(AlO4H4) ;
  
  double w_cl       = Flux(Cl) ;
  
  double w_cah2sio4 = Flux(CaH2SiO4) ;
  double w_cah3sio4 = Flux(CaH3SiO4) ;
  
  double w_naco3    = Flux(NaCO3) ;
  double w_nahco3   = Flux(NaHCO3) ;
  
  double w_caco3    = Flux(CaCO3) ;
  double w_cahco3   = Flux(CaHCO3) ;
  
  double w_caso4    = Flux(CaSO4) ;
  double w_cahso4   = Flux(CaHSO4) ;

  
  /* Fluxes as element: C, Ca, Si ... */
  /* Compounds type I */
  double w_ca_l    = w_ca + w_caoh + w_caoh2 ;
  double w_si_l    = w_h2sio4 + w_h3sio4 + w_h4sio4 ;
  double w_na_l    = w_na + w_naoh ;
  double w_k_l     = w_k + w_koh ;
  double w_c_l     = w_co2 + w_h2co3 + w_hco3 + w_co3 ;
  double w_s_l     = w_h2so4 + w_hso4 + w_so4 ;
  double w_al_l    = w_al + w_alo4h4 ;
  double w_cl_l    = w_cl ;
  /* Compounds type II */
  double w_ca_si_l = w_cah2sio4 + w_cah3sio4 ;
  double w_ca_c_l  = w_cahco3 + w_caco3 ;
  double w_ca_s_l  = w_caso4 + w_cahso4 ;
  double w_na_c_l  = w_nahco3 + w_naco3 ;
  
  
  /* Backup */
  
  CementSolutionDiffusion_GetElementFluxOf(csd,Ca) = w_ca_l + w_ca_si_l + w_ca_c_l + w_ca_s_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,Si) = w_si_l + w_ca_si_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,Na) = w_na_l + w_na_c_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,K)  = w_k_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,C)  = w_c_l + w_ca_c_l + w_na_c_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,S)  = w_s_l + w_ca_s_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,Al) = w_al_l ;
  CementSolutionDiffusion_GetElementFluxOf(csd,Cl) = w_cl_l ;
}


