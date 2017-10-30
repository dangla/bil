#ifndef ELEMENTSSOL_H
#define ELEMENTSSOL_H

/* class-like structures "ElementsSol_t" and attributes */

/* vacuous declarations and typedef names */
struct ElementsSol_s  ; typedef struct ElementsSol_s  ElementsSol_t ;
struct ElementSol_s   ; typedef struct ElementSol_s   ElementSol_t ;



/* Declaration of Macros, Methods and Structures */


/* 1. ElementsSol_t
 * -------------*/
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



/* 2. ElementSol_t 
 * ---------------*/
 
extern ElementSol_t* (ElementSol_GetDeepElementSol)(ElementSol_t*,unsigned int) ;


#define ElementSol_GetPreviousElementSol(ES)    ((ES)->prev)
#define ElementSol_GetImplicitGenericData(ES)   ((ES)->impgdat)
#define ElementSol_GetExplicitGenericData(ES)   ((ES)->expgdat)
#define ElementSol_GetConstantGenericData(ES)   ((ES)->cstgdat)




/* Nb of (im/ex)plicit and constant terms */
#define ElementSol_GetNbOfImplicitTerms(ES) \
        GenericData_GetNbOfData(ElementSol_GetImplicitGenericData(ES))

#define ElementSol_GetNbOfExplicitTerms(ES) \
        GenericData_GetNbOfData(ElementSol_GetExplicitGenericData(ES))

#define ElementSol_GetNbOfConstantTerms(ES) \
        GenericData_GetNbOfData(ElementSol_GetConstantGenericData(ES))




/* Access to (im/ex)plicit and constant terms */
#define ElementSol_GetImplicitTerm(ES) \
        GenericData_GetData(ElementSol_GetImplicitGenericData(ES))

#define ElementSol_GetExplicitTerm(ES) \
        GenericData_GetData(ElementSol_GetExplicitGenericData(ES))

#define ElementSol_GetConstantTerm(ES) \
        GenericData_GetData(ElementSol_GetConstantGenericData(ES))






struct ElementsSol_s {
  unsigned int NbOfElements ;
  unsigned int NbOfImpltTerms ;         /* nb of implicit terms */
  unsigned int NbOfExpltTerms ;         /* nb of explicit terms */
  unsigned int NbOfConstTerms ;         /* Nb of constant terms */
  ElementSol_t* elementsol ;
} ;


#include "GenericData.h"

struct ElementSol_s {              /* Element Solutions */
  GenericData_t* impgdat ;         /* Implicit generic data */
  GenericData_t* expgdat ;         /* Explicit generic data */
  GenericData_t* cstgdat ;         /* Constant generic data */
  ElementSol_t*  prev ;            /* Previous Element Solutions */
} ;


#endif
