#ifndef SHAPEFCTS_H
#define SHAPEFCTS_H


/* vacuous declarations and typedef names */

/* class-like structure "DataSet_t and attributes */
struct ShapeFcts_s      ; typedef struct ShapeFcts_s      ShapeFcts_t ;


extern ShapeFcts_t*  (ShapeFcts_Create)(void) ;
extern void          (ShapeFcts_Delete)(void*) ;
extern int           (ShapeFcts_FindShapeFct)(ShapeFcts_t*,int,int) ;
extern int           (ShapeFcts_AddShapeFct)(ShapeFcts_t*,int,int) ;


#define ShapeFcts_MaxNbOfShapeFcts             (10)

#define ShapeFcts_GetNbOfShapeFcts(SFS)    ((SFS)->n_sh)
#define ShapeFcts_GetShapeFct(SFS)         ((SFS)->sh)


#include "ShapeFct.h"

struct ShapeFcts_s {          /* Shape functions */
  unsigned int n_sh ;         /* Number of shape functions */
  ShapeFct_t*  sh ;           /* Shape function */
} ;





#endif
