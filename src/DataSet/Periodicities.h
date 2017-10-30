#ifndef PERIODICITIES_H
#define PERIODICITIES_H

/* vacuous declarations and typedef names */

/* class-like structures */
struct Periodicities_s       ; typedef struct Periodicities_s      Periodicities_t ;
/*   1. Periodicities_t attributes */
struct Periodicity_s         ; typedef struct Periodicity_s        Periodicity_t ;



/* 1. Periodicities_t 
 * ------------------*/
#include "DataFile.h"
#include "Mesh.h"
#include "Graph.h"

extern Periodicities_t* (Periodicities_Create)(DataFile_t*) ;
extern void             (Periodicities_ResetMatrixPermutationNumbering)(Mesh_t*) ;
extern void             (Periodicities_UpdateMatrixPermutationNumbering)(Mesh_t*) ;
extern void             (Periodicities_UpdateGraph)(Mesh_t*,Graph_t*) ;


#define Periodicities_GetNbOfPeriodicities(PS) ((PS)->nbperiod)
#define Periodicities_GetPeriodicity(PS)       ((PS)->periodicity)



/* 2. Periodicity_t 
 * ----------------*/
#define Periodicity_GetMasterRegion(P)         ((P)->masterreg)
#define Periodicity_GetSlaveRegion(P)          ((P)->slavereg)
#define Periodicity_GetPeriodVector(P)         ((P)->vector)


struct Periodicities_s {            /* Periodicities */
  unsigned int   nbperiod ;         /* Nb of periodicities */
  Periodicity_t* periodicity ;      /* Periodicity */
} ;

struct Periodicity_s {              /* Periodicity */
  int     masterreg ;               /* Master region index */
  int     slavereg ;                /* Slave region index */
  double* vector ;                  /* Period vector */
} ;



#endif
