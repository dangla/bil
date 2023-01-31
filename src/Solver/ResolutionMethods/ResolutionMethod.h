#ifndef RESOLUTIONMETHOD_H
#define RESOLUTIONMETHOD_H


enum ResolutionMethod_e {     /* Type of resolution method */
  ResolutionMethod_CROUT,     /* Crout method */
  ResolutionMethod_SuperLU,   /* SuperLU method */
  ResolutionMethod_MA38,      /* MA38 method */
  ResolutionMethod_SuperLUMT,  /* SuperLUMT method (multi-threaded) */
  ResolutionMethod_SuperLUDist  /* SuperLUDist method (distributed memory) */
} ;


/* class-like structures "ResolutionMethod_t" and attributes */

/* vacuous declarations and typedef names */
typedef enum ResolutionMethod_e  ResolutionMethod_e ;
struct ResolutionMethod_s   ; typedef struct ResolutionMethod_s  ResolutionMethod_t ;



#include "Options.h"

extern ResolutionMethod_t* (ResolutionMethod_Create)(Options_t*) ;
extern void                (ResolutionMethod_Delete)(void*) ;



/** The getters */
#define ResolutionMethod_GetType(RM)       ((RM)->type)
#define ResolutionMethod_GetOptions(RM)    ((RM)->options)


#include "Utils.h"

#define ResolutionMethod_Is(RM,KEY) \
        (ResolutionMethod_GetType(RM) == Utils_CAT(ResolutionMethod_,KEY))

#define ResolutionMethod_Type(KEY) \
        ((ResolutionMethod_e) Utils_CAT(ResolutionMethod_,KEY))



/* complete the structure types by using the typedef */

struct ResolutionMethod_s {
  ResolutionMethod_e type ;
  Options_t* options ;
} ;

#endif
