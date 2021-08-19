#ifndef ELEMENTSOL_H
#define ELEMENTSOL_H

/* class-like structures "ElementsSol_t" and attributes */

/* vacuous declarations and typedef names */
struct ElementSol_s   ; typedef struct ElementSol_s   ElementSol_t ;



 
extern ElementSol_t* (ElementSol_GetDeepElementSol)(ElementSol_t*,unsigned int) ;


#define ElementSol_GetPreviousElementSol(ES)    ((ES)->prev)
#define ElementSol_GetNextElementSol(ES)        ((ES)->next)
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


/* Add (im/ex)plicit and constant generic data */
#define ElementSol_AddImplicitGenericData(ES,GD) \
        GenericData_Append(ElementSol_GetImplicitGenericData(ES),GD)
        
#define ElementSol_AddExplicitGenericData(ES,GD) \
        GenericData_Append(ElementSol_GetExplicitGenericData(ES),GD)
        
#define ElementSol_AddConstantGenericData(ES,GD) \
        GenericData_Append(ElementSol_GetConstantGenericData(ES),GD)



/* Find (im/ex)plicit and constant generic data */
#define ElementSol_FindImplicitGenericData(ES,...) \
        GenericData_Find(ElementSol_GetImplicitGenericData(ES),__VA_ARGS__)
        
#define ElementSol_FindExplicitGenericData(ES,...) \
        GenericData_Find(ElementSol_GetExplicitGenericData(ES),__VA_ARGS__)
        
#define ElementSol_FindConstantGenericData(ES,...) \
        GenericData_Find(ElementSol_GetConstantGenericData(ES),__VA_ARGS__)



/* Find (im/ex)plicit and constant data */
#define ElementSol_FindImplicitData(ES,...) \
        GenericData_FindData(ElementSol_GetImplicitGenericData(ES),__VA_ARGS__)
        
#define ElementSol_FindExplicitData(ES,...) \
        GenericData_FindData(ElementSol_GetExplicitGenericData(ES),__VA_ARGS__)
        
#define ElementSol_FindConstantData(ES,...) \
        GenericData_FindData(ElementSol_GetConstantGenericData(ES),__VA_ARGS__)


#include "GenericData.h"


struct ElementSol_s {              /* Element Solutions */
  GenericData_t* impgdat ;         /* Implicit generic data */
  GenericData_t* expgdat ;         /* Explicit generic data */
  GenericData_t* cstgdat ;         /* Constant generic data */
  ElementSol_t*  prev ;            /* Previous Element Solutions */
  ElementSol_t*  next ;            /* Next Element Solutions */
} ;


#endif
