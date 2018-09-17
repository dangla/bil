#ifndef COORSYS_H
#define COORSYS_H

enum CoorSys_e {              /* coordinate system */
  COORSYS_CARTESIAN,
  COORSYS_CYLINDRICAL,
  COORSYS_SPHERICAL,
  CARTESIAN = COORSYS_CARTESIAN /* Should be eliminated */
} ;


/* vacuous declarations and typedef names */

/*     1. CoorSys_t attributes */
typedef enum CoorSys_e      CoorSys_t ;



/* Test the coordinate system */
#define CoorSys_IsCartesian(COO) \
        (COO == COORSYS_CARTESIAN)
        
#define CoorSys_IsCylindrical(COO) \
        (COO == COORSYS_CYLINDRICAL)
        
#define CoorSys_IsSpherical(COO) \
        (COO == COORSYS_SPHERICAL)



/* Set the coordinate system */
#define CoorSys_SetCartesian(COO) \
        (COO = COORSYS_CARTESIAN)
        
#define CoorSys_SetCylindrical(COO) \
        (COO = COORSYS_CYLINDRICAL)
        
#define CoorSys_SetSpherical(COO) \
        (COO = COORSYS_SPHERICAL)
        

#endif
