#ifndef ICOND_H
#define ICOND_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct ICond_s        ; typedef struct ICond_s        ICond_t ;


#define ICond_MaxLengthOfKeyWord        (30)
#define ICond_MaxLengthOfFileName       (60)

#define ICond_GetRegionIndex(ICD)             ((ICD)->reg)
#define ICond_GetNameOfUnknown(ICD)           ((ICD)->inc)
#define ICond_GetFunction(ICD)                ((ICD)->fn)
#define ICond_GetField(ICD)                   ((ICD)->ch)
#define ICond_GetFileNameOfNodalValues(ICD)   ((ICD)->file)



#include "Function.h"
#include "Field.h"


struct ICond_s {              /* Initial condition */
  int    reg ;                /* Region index */
  char*  inc ;                /* Name of unknown */
  Function_t* fn ;            /* Time function */
  Field_t* ch ;               /* Field */
  char*  file ;               /* Name of file containing nodal values */
} ;

#endif
