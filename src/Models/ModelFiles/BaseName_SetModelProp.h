#ifndef BASENAME_SETMODELPROP_H
#define BASENAME_SETMODELPROP_H


#include "Utils.h"

#define BASENAME_SETMODELPROP(base)   (Utils_CAT(base,_SetModelProp))

/* The macro BASENAME is sent from the compiler */
#if defined(BASENAME)
  #define BaseName_SetModelProp         BASENAME_SETMODELPROP(BASENAME)
#elif defined(__BASEFILE__)
  #define BaseName_SetModelProp         BASENAME_SETMODELPROP(__BASEFILE__)
#else
  #error "The macro BASENAME is not defined!"
#endif

#endif
