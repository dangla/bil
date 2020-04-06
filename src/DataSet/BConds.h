#ifndef BCONDS_H
#define BCONDS_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct BConds_s       ; typedef struct BConds_s       BConds_t ;


#include "DataFile.h"
#include "Functions.h"
#include "Fields.h"
#include "Mesh.h"


extern BConds_t* BConds_New(const int) ;
extern BConds_t* BConds_Create(DataFile_t*,Fields_t*,Functions_t*) ;
extern void      BConds_EliminateMatrixRowColumnIndexes(BConds_t*,Mesh_t*) ;
extern void      BConds_AssignBoundaryConditions(BConds_t*,Mesh_t*,double) ;


#define BConds_GetNbOfBConds(BCS)        ((BCS)->n_cl)
#define BConds_GetBCond(BCS)             ((BCS)->cl)


#include "BCond.h"


struct BConds_s {             /* boundary conditions */
  unsigned int n_cl ;         /* nb */
  BCond_t* cl ;               /* boundary condition */
} ;

#endif
