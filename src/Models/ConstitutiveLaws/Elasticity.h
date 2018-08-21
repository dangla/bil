#ifndef ELASTICITY_H
#define ELASTICITY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Elasticity_s     ; typedef struct Elasticity_s     Elasticity_t ;


/* 1. Elasticity_t */

extern Elasticity_t*  (Elasticity_Create)(void) ;
extern void           (Elasticity_Delete)(void*) ;


/* Accessors */
#define Elasticity_GetStiffnessTensor(EL)      ((EL)->cijkl)



struct Elasticity_s {
  double* cijkl ;
} ;

#endif
