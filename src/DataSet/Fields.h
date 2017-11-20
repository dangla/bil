#ifndef FIELDS_H
#define FIELDS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Fields_s       ; typedef struct Fields_s       Fields_t ;


#include "Materials.h"
#include "DataFile.h"
#include "Geometry.h"


extern Fields_t* Fields_Create(DataFile_t*,Materials_t*,Geometry_t*) ;


#define Fields_GetNbOfFields(FLDS)    ((FLDS)->n_ch)
#define Fields_GetField(FLDS)         ((FLDS)->ch)


#include "Field.h"


struct Fields_s {             /* fields */
  unsigned int n_ch ;         /* nb of fields */
  Field_t* ch ;               /* field */
} ;

#endif
