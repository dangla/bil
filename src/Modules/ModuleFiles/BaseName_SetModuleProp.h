#ifndef BASENAME_SETMODULEPROP_H
#define BASENAME_SETMODULEPROP_H


#include "Utils.h"

#define BASENAME_SETMODULEPROP(base)   (Utils_CAT(base,_SetModuleProp))

/* The macro BASENAME is sent from the compiler */
#define BaseName_SetModuleProp         BASENAME_SETMODULEPROP(BASENAME)

#endif
