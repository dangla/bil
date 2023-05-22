#ifndef ENTRY_H
#define ENTRY_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Entry_s  ; typedef struct Entry_s  Entry_t ;



extern Entry_t*    (Entry_Create)   (int,char**) ;
extern void        (Entry_Delete)   (void*) ;
extern int         (Entry_Execute)  (Entry_t*) ;


#define Entry_GetContext(E)          ((E)->context)



#include "Context.h"

struct Entry_s {
  Context_t*     context ;
} ;


#endif
