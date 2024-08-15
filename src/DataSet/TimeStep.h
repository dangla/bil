#ifndef TIMESTEP_H
#define TIMESTEP_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct TimeStep_s     ; typedef struct TimeStep_s     TimeStep_t ;


#include "DataFile.h"
#include "ObVals.h"
#include "Nodes.h"
#include "Solution.h"

extern TimeStep_t*  (TimeStep_Create)(DataFile_t*,ObVals_t*) ;
extern void         (TimeStep_Delete)(void*) ;
extern double       (TimeStep_ComputeTimeStep)(TimeStep_t*,Solution_t*,double,double) ;


#define TimeStep_GetInitialTimeStep(TS)       ((TS)->dtini)
#define TimeStep_GetMaximumTimeStep(TS)       ((TS)->dtmax)
#define TimeStep_GetMinimumTimeStep(TS)       ((TS)->dtmin)
#define TimeStep_GetMaximumCommonRatio(TS)    ((TS)->raison)
#define TimeStep_GetReductionFactor(TS)       ((TS)->fr)
#define TimeStep_GetObVals(TS)                ((TS)->obvals)
#define TimeStep_GetLocation(TS)              ((TS)->loc)
#define TimeStep_GetDataFile(TS)              ((TS)->datafile)
//#define TimeStep_GetSequentialIndex(TS)       ((TS)->sequentialindex)



/* Accesss to Objective Variations */
#define TimeStep_GetObVal(TS) \
        ObVals_GetObVal(TimeStep_GetObVals(TS))


/* Time location management */
#define TimeStep_SetLocationAtBegin(TS) \
        do {TimeStep_GetLocation(TS) = 0 ;} while(0)
        
#define TimeStep_SetLocationInBetween(TS) \
        do {TimeStep_GetLocation(TS) = 1 ;} while(0)
        
#define TimeStep_SetLocationAtEnd(TS) \
        do {TimeStep_GetLocation(TS) = 2 ;} while(0)
        
#define TimeStep_IsLocatedAtBegin(TS) \
        (TimeStep_GetLocation(TS) == 0)
        
#define TimeStep_IsLocatedInBetween(TS) \
        (TimeStep_GetLocation(TS) == 1)
        
#define TimeStep_IsLocatedAtEnd(TS) \
        (TimeStep_GetLocation(TS) == 2)
        
        

/* The sequential index */
#define TimeStep_GetSequentialIndex(TS) \
        DataFile_GetSequentialIndex(TimeStep_GetDataFile(TS))



struct TimeStep_s {           /* Time step management */
  double dtini ;              /* Initial time step */
  double dtmax ;              /* Maximum time step */
  double dtmin ;              /* Minimum time step */
  double raison ;             /* Maximum common ratio */
  double fr ;                 /* Factor reducing the time step */
  ObVals_t* obvals ;          /* Objective variations */
  char   loc ;                /* Time location */
  DataFile_t*    datafile ;   /* data file */
  //int sequentialindex ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
