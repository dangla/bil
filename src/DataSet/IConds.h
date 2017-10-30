#ifndef ICONDS_H
#define ICONDS_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct IConds_s       ; typedef struct IConds_s       IConds_t ;
/*   1. IConds_t attributes */
struct ICond_s        ; typedef struct ICond_s        ICond_t ;


/* 1. IConds_t
 * -----------*/
#include "DataFile.h"
#include "Functions.h"
#include "Mesh.h"
#include "Fields.h"

//extern IConds_t* IConds_Create(DataFile_t*,Mesh_t*,Fields_t*) ;
extern IConds_t* IConds_Create(DataFile_t*,Fields_t*,Functions_t*) ;
extern void      IConds_AssignInitialConditions(IConds_t*,Mesh_t*,double) ;

#define IConds_MaxLengthOfKeyWord        (30)
#define IConds_MaxLengthOfFileName       (60)


#define IConds_GetNbOfIConds(ICDS)              ((ICDS)->n_ic)
#define IConds_GetICond(ICDS)                   ((ICDS)->ic)
#define IConds_GetFileNameOfNodalValues(ICDS)   ((ICDS)->file)



/* 2. ICond_t
 * ----------*/
#define ICond_MaxLengthOfKeyWord        (30)
#define ICond_MaxLengthOfFileName       (60)

#define ICond_GetRegionIndex(ICD)             ((ICD)->reg)
#define ICond_GetNameOfUnknown(ICD)           ((ICD)->inc)
#define ICond_GetFunction(ICD)                ((ICD)->fn)
#define ICond_GetField(ICD)                   ((ICD)->ch)
#define ICond_GetFileNameOfNodalValues(ICD)   ((ICD)->file)



struct IConds_s {             /* Initial Conditions */
  char*  file ;               /* Name of file containing nodal values */
  unsigned int n_ic ;         /* Nb of IC */
  ICond_t* ic ;               /* Initial condition */
} ;


struct ICond_s {              /* Initial condition */
  int    reg ;                /* Region index */
  char*  inc ;                /* Name of unknown */
  Function_t* fn ;            /* Time function */
  Field_t* ch ;               /* Field */
  char*  file ;               /* Name of file containing nodal values */
} ;

#endif
