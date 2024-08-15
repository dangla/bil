#ifndef ELEMENTS_H
#define ELEMENTS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Elements_s     ; typedef struct Elements_s     Elements_t ;


#include "Materials.h"

extern Elements_t*  (Elements_New)                      (const int,const int) ;
extern void         (Elements_Delete)                   (void*) ;
extern void         (Elements_CreateMore)               (Elements_t*) ;
extern void         (Elements_LinkUp)                   (Elements_t*,Materials_t*) ;
extern void         (Elements_DefineProperties)         (Elements_t*) ;
extern int          (Elements_ComputeNbOfMatrixEntries) (Elements_t*) ;
extern int          (Elements_ComputeNbOfSelectedMatrixEntries)(Elements_t*,const int) ;
extern void         (Elements_UpdateMatrixRowColumnIndexes)(Elements_t*) ;
extern void         (Elements_EliminateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t*) ;
extern void         (Elements_UpdateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t*) ;
extern void         (Elements_InitializeMatrixRowColumnIndexes)(Elements_t*) ;


/* Accessors */
#define Elements_GetNbOfElements(ELTS)           ((ELTS)->NbOfElements)
#define Elements_GetNbOfConnectivities(ELTS)     ((ELTS)->NbOfConnectivities)
#define Elements_GetElement(ELTS)                ((ELTS)->Element)
#define Elements_GetPointerToNode(ELTS)          ((ELTS)->PointerToNode)
#define Elements_GetRegions(ELTS)                ((ELTS)->Regions)
#define Elements_GetShapeFcts(ELTS)              ((ELTS)->ShapeFcts)
#define Elements_GetIntFcts(ELTS)                ((ELTS)->IntFcts)
#define Elements_GetBuffers(ELTS)                ((ELTS)->Buffers)
#define Elements_GetMaximumSizeOfElements(ELTS)  ((ELTS)->MaximumSizeOfElements)
#define Elements_GetMinimumSizeOfElements(ELTS)  ((ELTS)->MinimumSizeOfElements)




#include "IntFcts.h"
#include "ShapeFcts.h"
#include "Buffers.h"
#include "Element.h"
#include "Node.h"
#include "Regions.h"


struct Elements_s {
  Element_t* Element ;
  Node_tt* PointerToNode ;
  Regions_t*   Regions ;
  ShapeFcts_t* ShapeFcts ;    /* Shape functions */
  IntFcts_t*   IntFcts ;      /* Interpolation functions */
  Buffers_t*  Buffers ;
  double      MaximumSizeOfElements ;
  double      MinimumSizeOfElements ;
  unsigned int NbOfElements ;
  unsigned int NbOfConnectivities ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
