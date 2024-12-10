#ifndef GEOMETRY_H
#define GEOMETRY_H

#ifdef __CPLUSPLUS
extern "C" {
#endif


/* vacuous declarations and typedef names */

/* class-like structure */
struct Geometry_s     ; typedef struct Geometry_s     Geometry_t ;


/* 1. Geometry_t 
 * -------------*/
#include "DataFile.h"

extern Geometry_t*  (Geometry_Create)(DataFile_t*) ;
extern void         (Geometry_Delete)(void*) ;
extern Geometry_t*  (Geometry_New)   (void) ;


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
          Symmetry_SetNoSymmetry(Geometry_GetSymmetry(GEO)) ; \
          CoorSys_SetCartesian(Geometry_GetCoordinateSystem(GEO)) ; \
        } while(0)
        
#define Geometry_SetPlaneSymmetry(GEO) \
        do { \
          Symmetry_SetPlaneSymmetry(Geometry_GetSymmetry(GEO)) ; \
          CoorSys_SetCartesian(Geometry_GetCoordinateSystem(GEO)) ; \
        } while(0)

#define Geometry_SetCylindricalSymmetry(GEO) \
        do { \
          Symmetry_SetCylindricalSymmetry(Geometry_GetSymmetry(GEO)) ; \
          CoorSys_SetCylindrical(Geometry_GetCoordinateSystem(GEO)) ; \
        } while(0)

#define Geometry_SetSphericalSymmetry(GEO) \
        do { \
          Symmetry_SetSphericalSymmetry(Geometry_GetSymmetry(GEO)) ; \
          CoorSys_SetSpherical(Geometry_GetCoordinateSystem(GEO)) ; \
        } while(0)


#include "Symmetry.h"
#include "CoorSys.h"
#include "Periodicities.h"

struct Geometry_s {
  unsigned int dim ;    /* Dimension (1,2,3) */
  Symmetry_t symmetry ;       /* Symmetry */
  CoorSys_t coorsys ;         /* Coordinate system */
  Periodicities_t* periodicities ;
} ;



#ifdef __CPLUSPLUS
}
#endif
#endif
