#ifndef UNIT_H
#define UNIT_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Unit_s         ; typedef struct Unit_s         Unit_t ;



#include "DataFile.h"

extern Unit_t* (Unit_New)     (void) ;
extern void    (Unit_Delete)  (void*) ;
extern int     (Unit_Scan)    (Unit_t*,DataFile_t*) ;



#define Unit_MaxLengthOfKeyWord   (30)



#define Unit_GetName(U)         ((U)->name)
#define Unit_GetValue(U)        ((U)->value)


struct Unit_s {           /* unit */
  char*  name ;           /* Its name */
  double value ;          /* Its value */
} ;


#endif
