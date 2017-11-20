#ifndef ICONDS_H
#define ICONDS_H

/* vacuous declarations and typedef names */

/* class-like structure and */
struct IConds_s       ; typedef struct IConds_s       IConds_t ;


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



#include "ICond.h"


struct IConds_s {             /* Initial Conditions */
  char*  file ;               /* Name of file containing nodal values */
  unsigned int n_ic ;         /* Nb of IC */
  ICond_t* ic ;               /* Initial condition */
} ;


#endif
