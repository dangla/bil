#ifndef FIELDS_H
#define FIELDS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct Fields_s       ; typedef struct Fields_s       Fields_t ;


//#include "Materials.h"
#include "Field.h"
#include "DataFile.h"


extern Fields_t*   (Fields_New)     (const int) ;
extern Fields_t*   (Fields_Create)  (DataFile_t*) ;
extern void        (Fields_Delete)  (void*) ;


#define Fields_GetNbOfFields(FLDS)    ((FLDS)->n_ch)
#define Fields_GetField(FLDS)         ((FLDS)->ch)




struct Fields_s {             /* fields */
  unsigned int n_ch ;         /* nb of fields */
  Field_t* ch ;               /* field */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
