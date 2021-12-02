#ifndef PERIODICITY_H
#define PERIODICITY_H

/* vacuous declarations and typedef names */

/* class-like structures */
struct Periodicity_s         ; typedef struct Periodicity_s        Periodicity_t ;


extern Periodicity_t* (Periodicity_New)(void) ;
extern void           (Periodicity_Delete)(void*) ;


#define Periodicity_GetMasterRegion(P)         ((P)->masterreg)
#define Periodicity_GetSlaveRegion(P)          ((P)->slavereg)
#define Periodicity_GetPeriodVector(P)         ((P)->vector)



struct Periodicity_s {              /* Periodicity */
  int     masterreg ;               /* Master region index */
  int     slavereg ;                /* Slave region index */
  double* vector ;                  /* Period vector */
} ;



#endif
