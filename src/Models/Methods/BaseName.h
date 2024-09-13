#ifndef BASENAME_H
#define BASENAME_H


#include "Utils.h"


/* The macro BASENAME is sent from the compiler */
#if defined(BASENAME)
  #define BaseName(method)     Utils_CAT(BASENAME,method)
#elif defined(__BASEFILE__)
  #define BaseName(method)     Utils_CAT(__BASEFILE__,method)
#else
  #error "The macros BASENAME or __BASEFILE__ are not defined!"
#endif

#endif
