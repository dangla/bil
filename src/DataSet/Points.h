#ifndef POINTS_H
#define POINTS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Points_s       ; typedef struct Points_s       Points_t ;


#include "DataFile.h"
#include "Mesh.h"


extern Points_t*  (Points_New)     (const int) ;
extern Points_t*  (Points_Create)  (DataFile_t*,Mesh_t*) ;
extern void       (Points_Delete)  (void*) ;



#define Points_GetNbOfPoints(PTS)    ((PTS)->n_points)
#define Points_GetPoint(PTS)         ((PTS)->point)



#include "Point.h"

struct Points_s {
  int n_points ;              /* nb of points */
  Point_t*  point ;           /* Point */
} ;


#endif
