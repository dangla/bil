#ifndef BASENAME_SETMODELPROP_H
#define BASENAME_SETMODELPROP_H

#define XBASENAME_SETMODELPROP(base)  (base##_SetModelProp)
#define BASENAME_SETMODELPROP(base)   XBASENAME_SETMODELPROP(base)

/* The macro BASENAME is sent from the compiler */
#define BaseName_SetModelProp         BASENAME_SETMODELPROP(BASENAME)

#endif
