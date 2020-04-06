#ifndef ELEMENTS_H
#define ELEMENTS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Elements_s     ; typedef struct Elements_s     Elements_t ;


#include "Materials.h"

extern Elements_t*  (Elements_Create)                   (const int) ;
extern Elements_t*  (Elements_New)                      (const int,const int) ;
extern void         (Elements_CreateMore)               (Elements_t*) ;
extern void         (Elements_LinkUp)                   (Elements_t*,Materials_t*) ;
extern void         (Elements_DefineProperties)         (Elements_t*) ;
extern int          (Elements_ComputeNbOfMatrixEntries) (Elements_t*) ;
extern void         (Elements_UpdateMatrixRowColumnIndexes)(Elements_t*) ;
extern void         (Elements_EliminateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t*) ;
extern void         (Elements_UpdateMatrixRowColumnIndexesOfOverlappingNodes)(Elements_t*) ;


/* Accessors */
#define Elements_GetNbOfElements(ELTS)           ((ELTS)->n_el)
#define Elements_GetNbOfConnectivities(ELTS)     ((ELTS)->n_con)
#define Elements_GetElement(ELTS)                ((ELTS)->el)
#define Elements_GetShapeFcts(ELTS)              ((ELTS)->shapefcts)
#define Elements_GetIntFcts(ELTS)                ((ELTS)->intfcts)
#define Elements_GetMaximumSizeOfElements(ELTS)  ((ELTS)->hmax)
#define Elements_GetMinimumSizeOfElements(ELTS)  ((ELTS)->hmin)





#include "IntFcts.h"
#include "ShapeFcts.h"
#include "Element.h"


struct Elements_s {           /* Elements */
  unsigned int n_el ;         /* Nb of elements */
  unsigned int n_con ;        /* Nb of connectivities */
  Element_t* el ;             /* Element */
  ShapeFcts_t* shapefcts ;    /* Shape functions */
  IntFcts_t* intfcts ;        /* Interpolation functions */
  double      hmax ;
  double      hmin ;
} ;

#endif
