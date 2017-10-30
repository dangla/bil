#ifndef BASENAME_SETMODULEPROP_H
#define BASENAME_SETMODULEPROP_H

#define XBASENAME_SETMODULEPROP(base)  (base##_SetModuleProp)
#define BASENAME_SETMODULEPROP(base)   XBASENAME_SETMODULEPROP(base)

/* The macro BASENAME is sent from the compiler */
#define BaseName_SetModuleProp         BASENAME_SETMODULEPROP(BASENAME)

#endif
