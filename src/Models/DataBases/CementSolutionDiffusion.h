#ifndef CEMENTSOLUTIONDIFFUSION_H
#define CEMENTSOLUTIONDIFFUSION_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct CementSolutionDiffusion_s     ; 
typedef struct CementSolutionDiffusion_s     CementSolutionDiffusion_t ;



extern CementSolutionDiffusion_t* (CementSolutionDiffusion_Create)(void) ;
extern CementSolutionDiffusion_t* (CementSolutionDiffusion_GetInstance)(void) ;
extern void   (CementSolutionDiffusion_ComputeFluxes)(CementSolutionDiffusion_t*) ;



#include "Elements.h"
#include "IntFcts.h"

#define CementSolutionDiffusion_MaxNbOfPotentialVectors \
        (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))



/* Accessors (Getters) */
#define CementSolutionDiffusion_GetTemperature(csc) \
      ((csc)->temperature)

#define CementSolutionDiffusion_GetDiffusionCoefficient(csd) \
      ((csd)->diffusioncoefficient)
      
#define CementSolutionDiffusion_GetGradient(csd) \
      ((csd)->gradient)
      
#define CementSolutionDiffusion_GetPotential(csd) \
      ((csd)->potential)
      
#define CementSolutionDiffusion_GetPointerToPotentials(csd) \
      ((csd)->pointertopotentials)
      
#define CementSolutionDiffusion_GetFlux(csd) \
      ((csd)->flux)
      
#define CementSolutionDiffusion_GetElementFlux(csd) \
      ((csd)->elementflux)
      
#define CementSolutionDiffusion_GetIonCurrent(csd) \
      ((csd)->ioncurrent)



#include "CementSolutionChemistry.h"


/* Macros for diffusion coefficients
 * ---------------------------------*/
        
#define CementSolutionDiffusion_NbOfSpecies \
        CementSolutionChemistry_NbOfSpecies

#define CementSolutionDiffusion_GetDiffusionCoefficientOf(csd,CPD) \
       (CementSolutionDiffusion_GetDiffusionCoefficient(csd)[CementSolutionChemistry_GetIndexOf(CPD)])
       

/* Synomyms */
#define CementSolutionDiffusion_NbOfConcentrations \
        CementSolutionDiffusion_NbOfSpecies



/* Macros for the potentials
 * -------------------------*/
#define CementSolutionDiffusion_GetPotentialAtPoint(csd,i) \
       (CementSolutionDiffusion_GetPointerToPotentials(csd)[i])



/* Macros for the concentration fluxes
 * -----------------------------------*/
#define CementSolutionDiffusion_GetFluxOf(csd,CPD) \
       (CementSolutionDiffusion_GetFlux(csd)[CementSolutionChemistry_GetIndexOf(CPD)])



/* Macros for the element fluxes
 * -----------------------------*/
#define CementSolutionDiffusion_NbOfElementFluxes \
        CementSolutionChemistry_NbOfElementConcentrations

#define CementSolutionDiffusion_GetElementFluxOf(csd,A) \
       (CementSolutionDiffusion_GetElementFlux(csd)[CementSolutionChemistry_E_##A])
       



struct CementSolutionDiffusion_s {
  double  temperature ;
  double* diffusioncoefficient ;
  double* gradient ;
  double* potential ;
  double** pointertopotentials ;
  double* flux ;
  double* elementflux ;
  double  ioncurrent ;
} ;

#endif
