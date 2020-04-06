#ifndef BASENAME_SETMODELPROP_H
#define BASENAME_SETMODELPROP_H



#include "Utils.h"

#define BASENAME_SETMODELPROP(base)   (Utils_CAT(base,_SetModelProp))

/* The macro BASENAME is sent from the compiler */
#define BaseName_SetModelProp         BASENAME_SETMODELPROP(BASENAME)

#endif
