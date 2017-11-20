#ifndef ELEMENTSSOL_H
#define ELEMENTSSOL_H

/* class-like structures "ElementsSol_t" and attributes */

/* vacuous declarations and typedef names */
struct ElementsSol_s  ; typedef struct ElementsSol_s  ElementsSol_t ;



#include "Mesh.h"

extern ElementsSol_t*   (ElementsSol_Create)(Mesh_t*) ;
extern void             (ElementsSol_AllocateMemoryForImplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForExplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForConstantTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_ShareConstantTermsFrom)(ElementsSol_t*,ElementsSol_t*) ;
extern void             (ElementsSol_Copy)(ElementsSol_t*,ElementsSol_t*) ;


#define ElementsSol_GetNbOfImplicitTerms(ESS)    ((ESS)->NbOfImpltTerms)
#define ElementsSol_GetNbOfExplicitTerms(ESS)    ((ESS)->NbOfExpltTerms)
#define ElementsSol_GetNbOfConstantTerms(ESS)    ((ESS)->NbOfConstTerms)
#define ElementsSol_GetNbOfElements(ESS)         ((ESS)->NbOfElements)
#define ElementsSol_GetElementSol(ESS)           ((ESS)->elementsol)




#include "ElementSol.h"


struct ElementsSol_s {
  unsigned int NbOfElements ;
  unsigned int NbOfImpltTerms ;         /* nb of implicit terms */
  unsigned int NbOfExpltTerms ;         /* nb of explicit terms */
  unsigned int NbOfConstTerms ;         /* Nb of constant terms */
  ElementSol_t* elementsol ;
} ;


#endif
