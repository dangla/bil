#ifndef UNITS_H
#define UNITS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Units_s        ; typedef struct Units_s        Units_t ;



#include "DataFile.h"
 
extern Units_t* Units_Create(DataFile_t*) ;


#include "Unit.h"


struct Units_s {              /* units */
  int i ;
} ;

#endif
