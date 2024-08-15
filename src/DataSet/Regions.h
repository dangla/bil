#ifndef REGIONS_H
#define REGIONS_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Regions_s        ; typedef struct Regions_s        Regions_t ;


#include "Region.h"

extern Regions_t*  (Regions_New)            (void) ;
extern void        (Regions_Delete)         (void*) ;
extern int         (Regions_FindRegionIndex)(Regions_t*,const char*) ;
extern Region_t*   (Regions_FindRegion)     (Regions_t*,const char*) ;


        
#define Regions_MaxNbOfRegions  (20)


#define Regions_GetNbOfRegions(RS)          ((RS)->NbOfRegions)
#define Regions_GetRegion(RS)               ((RS)->Region)



struct Regions_s {
  Region_t*    Region ;
  unsigned int NbOfRegions ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
