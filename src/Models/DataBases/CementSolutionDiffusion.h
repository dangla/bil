#ifndef CEMENTSOLUTIONDIFFUSION_H
#define CEMENTSOLUTIONDIFFUSION_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct CementSolutionDiffusion_s     ; 
typedef struct CementSolutionDiffusion_s     CementSolutionDiffusion_t ;



extern CementSolutionDiffusion_t* (CementSolutionDiffusion_Create)(void) ;
extern void   (CementSolutionDiffusion_ComputeFluxes)(CementSolutionDiffusion_t*) ;



#include "Elements.h"
#include "IntFcts.h"

#define CementSolutionDiffusion_MaxNbOfPotentialVectors \
        (MAX(Element_MaxNbOfNodes,IntFct_MaxNbOfIntPoints))



/* Accessors (Getters) */
#define CementSolutionDiffusion_GetTemperature(CSD) \
      ((CSD)->temperature)

#define CementSolutionDiffusion_GetDiffusionCoefficient(CSD) \
      ((CSD)->diffusioncoefficient)
      
#define CementSolutionDiffusion_GetGradient(CSD) \
      ((CSD)->gradient)
      
#define CementSolutionDiffusion_GetPotential(CSD) \
      ((CSD)->potential)
      
#define CementSolutionDiffusion_GetPointerToPotentials(CSD) \
      ((CSD)->pointertopotentials)
      
#define CementSolutionDiffusion_GetFlux(CSD) \
      ((CSD)->flux)
      
#define CementSolutionDiffusion_GetElementFlux(CSD) \
      ((CSD)->elementflux)
      
#define CementSolutionDiffusion_GetIonCurrent(CSD) \
      ((CSD)->ioncurrent)




/* Macros for the room temperature
 * -------------------------------*/
#define CementSolutionDiffusion_GetRoomTemperature(CSD) \
        Temperature_GetRoomValue(CementSolutionDiffusion_GetTemperature(CSD))
        

#define CementSolutionDiffusion_SetRoomTemperature(CSD,T) \
        Temperature_SetRoomTemperature(CementSolutionDiffusion_GetTemperature(CSD),T)





#include "CementSolutionChemistry.h"


/* Macros for diffusion coefficients
 * ---------------------------------*/
        
#define CementSolutionDiffusion_NbOfSpecies \
        CementSolutionChemistry_NbOfSpecies

#define CementSolutionDiffusion_GetDiffusionCoefficientOf(CSD,CPD) \
        (CementSolutionDiffusion_GetDiffusionCoefficient(CSD)[CementSolutionChemistry_GetIndexOf(CPD)])
       

/* Synomyms */
#define CementSolutionDiffusion_NbOfConcentrations \
        CementSolutionDiffusion_NbOfSpecies



/* Macros for the potentials
 * -------------------------*/
#define CementSolutionDiffusion_GetPotentialAtPoint(CSD,i) \
        (CementSolutionDiffusion_GetPointerToPotentials(CSD)[i])



/* Macros for the concentration fluxes
 * -----------------------------------*/
#define CementSolutionDiffusion_GetFluxOf(CSD,CPD) \
        (CementSolutionDiffusion_GetFlux(CSD)[CementSolutionChemistry_GetIndexOf(CPD)])



/* Macros for the element fluxes
 * -----------------------------*/
#define CementSolutionDiffusion_NbOfElementFluxes \
        CementSolutionChemistry_NbOfElementConcentrations

#define CementSolutionDiffusion_GetElementFluxOf(CSD,A) \
        (CementSolutionDiffusion_GetElementFlux(CSD)[CementSolutionChemistry_E_##A])
       


#include "Temperature.h"

struct CementSolutionDiffusion_s {
  Temperature_t* temperature ;
  double* diffusioncoefficient ;
  double* gradient ;
  double* potential ;
  double** pointertopotentials ;
  double* flux ;
  double* elementflux ;
  double  ioncurrent ;
} ;

#endif
