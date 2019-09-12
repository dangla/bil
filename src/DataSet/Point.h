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


#define Point_GetCoordinate(PT)                    ((PT)->coor)
#define Point_GetEnclosingElement(PT)              ((PT)->elt)
#define Point_GetRegionIndex(PT)                   ((PT)->reg)


#include "Element.h"

struct Point_s {
  double* coor ;              /* Coordinates (x,y,z) */
  int reg ;                   /* Region index */
  Element_t* elt ;            /* Element inside which the point lies */
} ;


#endif
