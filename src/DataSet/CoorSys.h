#ifndef COORSYS_H
#define COORSYS_H

enum CoorSys_e {              /* coordinate system */
  CoorSys_Cartesian,
  CoorSys_Cylindrical,
  CoorSys_Spherical,
  CARTESIAN = CoorSys_Cartesian /* Should be eliminated */
} ;


/* vacuous declarations and typedef names */

/*     1. CoorSys_t attributes */
typedef enum CoorSys_e      CoorSys_t ;



/* Test the coordinate system */
#define CoorSys_IsCartesian(COO) \
        (COO == CoorSys_Cartesian)
        
#define CoorSys_IsCylindrical(COO) \
        (COO == CoorSys_Cylindrical)
        
#define CoorSys_IsSpherical(COO) \
        (COO == CoorSys_Spherical)



/* Set the coordinate system */
#define CoorSys_SetCartesian(COO) \
        (COO = CoorSys_Cartesian)
        
#define CoorSys_SetCylindrical(COO) \
        (COO = CoorSys_Cylindrical)
        
#define CoorSys_SetSpherical(COO) \
        (COO = CoorSys_Spherical)
        

#endif
