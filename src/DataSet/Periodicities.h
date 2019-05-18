#ifndef PERIODICITIES_H
#define PERIODICITIES_H

/* vacuous declarations and typedef names */

/* class-like structures */
struct Periodicities_s       ; typedef struct Periodicities_s      Periodicities_t ;



#include "DataFile.h"
#include "Mesh.h"
#include "Graph.h"


extern Periodicities_t* (Periodicities_Create)(DataFile_t*) ;
extern Periodicities_t* (Periodicities_New)(const int) ;
extern void             (Periodicities_Delete)(void*) ;
extern void             (Periodicities_ResetMatrixPermutationNumbering)(Mesh_t*) ;
extern void             (Periodicities_UpdateMatrixPermutationNumbering)(Mesh_t*) ;
extern void             (Periodicities_UpdateGraph)(Mesh_t*,Graph_t*) ;


#define Periodicities_GetNbOfPeriodicities(PS) ((PS)->nbperiod)
#define Periodicities_GetPeriodicity(PS)       ((PS)->periodicity)



#include "Periodicity.h"


struct Periodicities_s {            /* Periodicities */
  unsigned int   nbperiod ;         /* Nb of periodicities */
  Periodicity_t* periodicity ;      /* Periodicity */
} ;


#endif
