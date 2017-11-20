#ifndef POINT_H
#define POINT_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Point_s        ; typedef struct Point_s        Point_t ;




#define Point_GetCoordinate(PT)                    ((PT)->coor)
#define Point_GetEnclosingElement(PT)              ((PT)->elt)


#include "Element.h"

struct Point_s {
  double* coor ;              /* Coordinates (x,y,z) */
  Element_t* elt ;            /* Element inside which the point lies */
} ;


#endif
