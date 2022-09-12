#ifndef BIL_H
#define BIL_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Bil_s  ; typedef struct Bil_s  Bil_t ;



extern Bil_t*    (Bil_Create)   (int,char**) ;
extern void      (Bil_Delete)   (void*) ;
extern int       (Bil_Execute)  (Bil_t*) ;


#define Bil_GetContext(BIL)          ((BIL)->context)



#include "Context.h"

struct Bil_s {
  Context_t*     context ;
} ;


#endif
