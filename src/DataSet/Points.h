#ifndef POINTS_H
#define POINTS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Points_s       ; typedef struct Points_s       Points_t ;
/*     1. Points_t attributes */
struct Point_s        ; typedef struct Point_s        Point_t ;


/* 1. Points_t 
 * -----------*/
#include "DataFile.h"
#include "Mesh.h"

extern Points_t*  Points_Create(DataFile_t*,Mesh_t*) ;

#define Points_GetNbOfPoints(PTS)    ((PTS)->n_points)
#define Points_GetPoint(PTS)         ((PTS)->point)



/* 2. Point_t 
 * ----------*/
#include "Elements.h"

#define Point_GetCoordinate(PT)                    ((PT)->coor)
#define Point_GetEnclosingElement(PT)              ((PT)->elt)


struct Points_s {
  unsigned int n_points ;     /* nb of points */
  Point_t*  point ;           /* Point */
} ;

struct Point_s {
  double* coor ;              /* Coordinates (x,y,z) */
  Element_t* elt ;            /* Element inside which the point lies */
} ;


#endif
