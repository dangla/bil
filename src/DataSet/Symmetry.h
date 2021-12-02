#ifndef SYMMETRY_H
#define SYMMETRY_H

enum Symmetry_e {             /* symmetry of the problem */
  Symmetry_None,
  Symmetry_Plane,
  Symmetry_Cylindrical,
  Symmetry_Spherical,
  AXIS = Symmetry_Cylindrical,       /* Should be eliminated */
  SPHE = Symmetry_Spherical, 
} ;


/* vacuous declarations and typedef names */

/*     1. Symmetry_t attributes */
typedef enum Symmetry_e     Symmetry_t ;



/* Test the symmetry */
#define Symmetry_IsCylindrical(SYM) \
        (SYM == Symmetry_Cylindrical)
        
#define Symmetry_IsSpherical(SYM) \
        (SYM == Symmetry_Spherical)
        
#define Symmetry_IsPlane(SYM) \
        (SYM == Symmetry_Plane)



/* Set the symmetry */
#define Symmetry_SetNoSymmetry(SYM) \
        (SYM = Symmetry_None)
        
#define Symmetry_SetPlaneSymmetry(SYM) \
        (SYM = Symmetry_Plane)
        
#define Symmetry_SetCylindricalSymmetry(SYM) \
        (SYM = Symmetry_Cylindrical)
        
#define Symmetry_SetSphericalSymmetry(SYM) \
        (SYM = Symmetry_Spherical)



/* Old notations which I try to eliminate little by little */
#define geom_t    Symmetry_t

#endif
