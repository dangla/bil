#ifndef ELEMENTS_H
#define ELEMENTS_H


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
#define Elements_GetNbOfElements(ELTS)           ((ELTS)->n_el)
#define Elements_GetNbOfConnectivities(ELTS)     ((ELTS)->n_con)
#define Elements_GetElement(ELTS)                ((ELTS)->el)
#define Elements_GetPointerToNode(ELTS)          ((ELTS)->node)
#define Elements_GetShapeFcts(ELTS)              ((ELTS)->shapefcts)
#define Elements_GetIntFcts(ELTS)                ((ELTS)->intfcts)
#define Elements_GetBuffers(ELTS)                ((ELTS)->buffers)
#define Elements_GetMaximumSizeOfElements(ELTS)  ((ELTS)->hmax)
#define Elements_GetMinimumSizeOfElements(ELTS)  ((ELTS)->hmin)





#include "IntFcts.h"
#include "ShapeFcts.h"
#include "Buffers.h"
#include "Element.h"
#include "Node.h"


struct Elements_s {           /* Elements */
  Element_t* el ;             /* Element */
  Node_tt* node ;             /* pointers to nodes */
  ShapeFcts_t* shapefcts ;    /* Shape functions */
  IntFcts_t* intfcts ;        /* Interpolation functions */
  Buffers_t*  buffers ;       /* Buffers */
  double      hmax ;
  double      hmin ;
  unsigned int n_el ;         /* Nb of elements */
  unsigned int n_con ;        /* Nb of connectivities */
} ;

#endif
