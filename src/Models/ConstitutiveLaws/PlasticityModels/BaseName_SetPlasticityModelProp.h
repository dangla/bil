#ifndef BASENAME_SETPLASTICITYMODELPROP_H
#define BASENAME_SETPLASTICITYMODELPROP_H


#include "Utils.h"

#define BASENAME_SETPLASTICITYMODELPROP(base)   (Utils_CAT(base,_SetModelProp))

/* The macro BASENAME is sent from the compiler */
#ifdef  BASENAME
  #define BaseName_SetPlasticityModelProp         BASENAME_SETPLASTICITYMODELPROP(BASENAME)
#else
  #error "The macro BASENAME is not defined!"
#endif

#endif
