#ifndef GEOMETRY_H
#define GEOMETRY_H

enum Symmetry_e {             /* symmetry of the problem */
  NO,
  PLAN,
  AXIS,
  SPHE
} ;

enum CoorSys_e {              /* coordinate system */
  CARTESIAN,
  CYLINDRICAL,
  SPHERICAL
} ;


/* vacuous declarations and typedef names */

/* class-like structure */
struct Geometry_s     ; typedef struct Geometry_s     Geometry_t ;

/*     1. Geometry_t attributes */
typedef enum CoorSys_e      CoorSys_t ;
typedef enum Symmetry_e     Symmetry_t ;



/* Test the symmetry */
#define Symmetry_IsCylindrical(SYM)     (SYM == AXIS)
#define Symmetry_IsSpherical(SYM)       (SYM == SPHE)
#define Symmetry_IsPlane(SYM)           (SYM == PLAN)



/* 1. Geometry_t 
 * -------------*/
#include "DataFile.h"

extern Geometry_t*  Geometry_Create(DataFile_t*) ;

#define Geometry_GetDimension(GEO)              ((GEO)->dim)
#define Geometry_GetSymmetry(GEO)               ((GEO)->symmetry)
#define Geometry_GetCoordinateSystem(GEO)       ((GEO)->coorsys)
#define Geometry_GetPeriodicities(GEO)          ((GEO)->periodicities)


/* Is it periodic? */
#define Geometry_IsPeriodic(GEO) \
        (Geometry_GetPeriodicities(GEO) != NULL)


/* Test the symmetry */
#define Geometry_HasCylindricalSymmetry(GEO) \
        Symmetry_IsCylindrical(Geometry_GetSymmetry(GEO))
        
#define Geometry_HasSphericalSymmetry(GEO) \
        Symmetry_IsSpherical(Geometry_GetSymmetry(GEO))


/* Set the symmetry */
#define Geometry_SetNoSymmetry(GEO) \
        do { \
          Geometry_GetSymmetry(GEO) = NO ; \
          Geometry_GetCoordinateSystem(GEO) = CARTESIAN ; \
        } while(0)
        
#define Geometry_SetPlaneSymmetry(GEO) \
        do { \
          Geometry_GetSymmetry(GEO) = PLAN ; \
          Geometry_GetCoordinateSystem(GEO) = CARTESIAN ; \
        } while(0)

#define Geometry_SetCylindricalSymmetry(GEO) \
        do { \
          Geometry_GetSymmetry(GEO) = AXIS ; \
          Geometry_GetCoordinateSystem(GEO) = CYLINDRICAL ; \
        } while(0)

#define Geometry_SetSphericalSymmetry(GEO) \
        do { \
          Geometry_GetSymmetry(GEO) = SPHE ; \
          Geometry_GetCoordinateSystem(GEO) = SPHERICAL ; \
        } while(0)


#include "Periodicities.h"

struct Geometry_s {
  unsigned short int dim ;    /* dimension (1,2,3) */
  Symmetry_t symmetry ;       /* symmetry (PLAN,AXIS,SPHE) */
  CoorSys_t coorsys ;         /* coordinate system (CARTESIAN,CYLINDRICAL,SPHERICAL) */
  Periodicities_t* periodicities ;
} ;


/* Old notations which I try to eliminate little by little */
#define geom_t    Symmetry_t

#endif
