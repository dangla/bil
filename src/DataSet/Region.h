#ifndef REGION_H
#define REGION_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Region_s         ; typedef struct Region_s         Region_t ;


extern Region_t*  (Region_New)   (void) ;
extern void       (Region_Delete)(void*) ;


        
#define Region_MaxLengthOfRegionName  (50)


/* Accessors */
//#define Region_GetRegionIndex(R)          ((R)->RegionIndex)
//#define Region_GetRegionTag(R)            ((R)->RegionTag)
#define Region_GetRegionName(R)           ((R)->RegionName)


struct Region_s {
  //int   RegionIndex ;
  //int   RegionTag ;
  char* RegionName ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
