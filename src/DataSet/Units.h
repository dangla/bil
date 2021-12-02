#ifndef UNITS_H
#define UNITS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Units_s        ; typedef struct Units_s        Units_t ;



#include "DataFile.h"
 
 
extern Units_t* (Units_New)     (void) ;
extern Units_t* (Units_Create)  (DataFile_t*) ;
extern void     (Units_Delete)  (void*) ;


#define Units_MaxNbOfUnits   (7)


#define Units_GetNbOfUnits(U)         ((U)->n_units)
#define Units_GetUnit(U)              ((U)->unit)


#include "Unit.h"


struct Units_s {              /* units */
  unsigned int n_units ;      /* nb of units */
  Unit_t*  unit ;             /* Unit */
} ;

#endif
