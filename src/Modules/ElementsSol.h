#ifndef ELEMENTSSOL_H
#define ELEMENTSSOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "ElementsSol_t" and attributes */

/* vacuous declarations and typedef names */
struct ElementsSol_s  ; typedef struct ElementsSol_s  ElementsSol_t ;



#include "Mesh.h"

extern ElementsSol_t*   (ElementsSol_Create)(Mesh_t*) ;
extern void             (ElementsSol_Delete)(void*) ;
extern void             (ElementsSol_AllocateMemoryForImplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForExplicitTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_AllocateMemoryForConstantTerms)(ElementsSol_t*) ;
extern void             (ElementsSol_Copy)(ElementsSol_t*,ElementsSol_t*) ;


//#define ElementsSol_GetNbOfImplicitTerms(ESS)    ((ESS)->NbOfImpltTerms)
//#define ElementsSol_GetNbOfExplicitTerms(ESS)    ((ESS)->NbOfExpltTerms)
//#define ElementsSol_GetNbOfConstantTerms(ESS)    ((ESS)->NbOfConstTerms)
#define ElementsSol_GetNbOfElements(ESS)         ((ESS)->NbOfElements)
#define ElementsSol_GetElementSol(ESS)           ((ESS)->elementsol)




/* Delete data */
#define ElementsSol_DeleteExplicitGenericData(ESS) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          int ElementsSol_i ; \
          for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_DeleteExplicitGenericData(elementsol + ElementsSol_i) ; \
          } \
        } while(0)
        
        
#define ElementsSol_DeleteConstantGenericData(ESS) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          int ElementsSol_i ; \
          for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_DeleteConstantGenericData(elementsol + ElementsSol_i) ; \
          } \
        } while(0)
        
        

/* Share data */
#define ElementsSol_ShareExplicitGenericData(ESS_DEST,ESS_SRC) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS_SRC) ; \
          ElementSol_t* elementsol_s = ElementsSol_GetElementSol(ESS_SRC) ; \
          ElementSol_t* elementsol_d = ElementsSol_GetElementSol(ESS_DEST) ; \
          { \
            int ElementsSol_i ; \
            for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
              GenericData_t* gdat = ElementSol_GetExplicitGenericData(elementsol_s + ElementsSol_i) ; \
              ElementSol_GetExplicitGenericData(elementsol_d + ElementsSol_i) = gdat ; \
            } \
          } \
        } while(0)


#define ElementsSol_ShareConstantGenericData(ESS_DEST,ESS_SRC) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS_SRC) ; \
          ElementSol_t* elementsol_s = ElementsSol_GetElementSol(ESS_SRC) ; \
          ElementSol_t* elementsol_d = ElementsSol_GetElementSol(ESS_DEST) ; \
          { \
            int ElementsSol_i ; \
            for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
              GenericData_t* gdat = ElementSol_GetConstantGenericData(elementsol_s + ElementsSol_i) ; \
              ElementSol_GetConstantGenericData(elementsol_d + ElementsSol_i) = gdat ; \
            } \
          } \
        } while(0)



/* Set pointers to null */
#define ElementsSol_DiscardExplicitGenericData(ESS) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          int ElementsSol_i ; \
          for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_GetExplicitGenericData(elementsol + ElementsSol_i) = NULL ; \
          } \
        } while(0)


#define ElementsSol_DiscardConstantGenericData(ESS) \
        do { \
          int NbOfElements = ElementsSol_GetNbOfElements(ESS) ; \
          ElementSol_t* elementsol = ElementsSol_GetElementSol(ESS) ; \
          int ElementsSol_i ; \
          for(ElementsSol_i = 0 ; ElementsSol_i < NbOfElements ; ElementsSol_i++) { \
            ElementSol_GetConstantGenericData(elementsol + ElementsSol_i) = NULL ; \
          } \
        } while(0)




#include "ElementSol.h"


struct ElementsSol_s {
  unsigned int NbOfElements ;
  //unsigned int NbOfImpltTerms ;         /* nb of implicit terms */
  //unsigned int NbOfExpltTerms ;         /* nb of explicit terms */
  //unsigned int NbOfConstTerms ;         /* Nb of constant terms */
  ElementSol_t* elementsol ;
} ;



#ifdef __CPLUSPLUS
}
#endif
#endif
