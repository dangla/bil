#ifndef UNITS_H
#define UNITS_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Units_s        ; typedef struct Units_s        Units_t ;
/*     1. Units_t attributes */
struct Unit_s         ; typedef struct Unit_s         Unit_t ;



/* Declaration of Macros, Methods and Structures */

/* 1. Units_t 
 * ----------*/

#include "DataFile.h"
 
extern Units_t* Units_Create(DataFile_t*) ;


#define Unit_MaxLengthOfKeyWord   (30)



struct Units_s {              /* units */
  int i ;
} ;

#endif
