#ifndef BCONDS_H
#define BCONDS_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct BConds_s       ; typedef struct BConds_s       BConds_t ;
/*   1. BConds_t attributes */
struct BCond_s        ; typedef struct BCond_s        BCond_t ;



/* 1. BConds_t
 * -----------*/
#include "DataFile.h"
#include "Functions.h"
#include "Fields.h"
#include "Mesh.h"

extern BConds_t* BConds_Create(DataFile_t*,Fields_t*,Functions_t*) ;
extern void      BConds_ResetMatrixNumbering(BConds_t*,Mesh_t*) ;
extern void      BConds_AssignBoundaryConditions(BConds_t*,Mesh_t*,double) ;

#define BConds_GetNbOfBConds(BCS)        ((BCS)->n_cl)
#define BConds_GetBCond(BCS)             ((BCS)->cl)



/* 2. BCond_t
 * ----------*/
#define BCond_MaxLengthOfKeyWord        (30)

#define BCond_GetRegionIndex(BC)         ((BC)->reg)
#define BCond_GetNameOfUnknown(BC)       ((BC)->inc)
#define BCond_GetNameOfEquation(BC)      ((BC)->eqn)
#define BCond_GetFunction(BC)            ((BC)->fn)
#define BCond_GetField(BC)               ((BC)->ch)

struct BConds_s {             /* boundary conditions */
  unsigned int n_cl ;         /* nb */
  BCond_t* cl ;               /* boundary condition */
} ;

struct BCond_s {              /* condition a la limite */
  int    reg ;                /* numero de la region */
  char*  inc ;                /* nom de l'inconnue a imposer */
  char*  eqn ;                /* nom de l'equation a eliminer */
  Function_t* fn ;            /* fonction du temps */
  Field_t* ch ;               /* champ */
} ;


/* Old notations which I try to eliminate little by little */
#define cond_t    BCond_t

#endif
