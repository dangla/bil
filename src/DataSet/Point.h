#ifndef POINT_H
#define POINT_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Point_s        ; typedef struct Point_s        Point_t ;


#include "Mesh.h"


extern Point_t*  (Point_New)                  (void) ;
extern void      (Point_Delete)               (void*) ;
extern void      (Point_SetEnclosingElement)  (Point_t*,Mesh_t*) ;
extern void      (Point_Scan)                 (Point_t*,char*) ;


#include "Region.h"

#define Point_MaxLengthOfRegionName      Region_MaxLengthOfRegionName


#define Point_GetCoordinate(PT)                    ((PT)->Coordinate)
#define Point_GetEnclosingElement(PT)              ((PT)->EnclosingElement)
#define Point_GetRegionTag(PT)                     ((PT)->RegionTag)
#define Point_GetRegionName(PT)                    ((PT)->RegionName)


#include "Element.h"

struct Point_s {
  double* Coordinate ;
  int   RegionTag ;
  char* RegionName ;
  Element_t* EnclosingElement ;  /* Element inside which the point lies */
} ;


#endif
