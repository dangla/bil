#ifndef SYMMETRY_H
#define SYMMETRY_H

enum Symmetry_e {             /* symmetry of the problem */
  SYMMETRY_NO,
  SYMMETRY_PLAN,
  SYMMETRY_AXIS,
  SYMMETRY_SPHE,
  AXIS = SYMMETRY_AXIS,       /* Should be eliminated */
  SPHE = SYMMETRY_SPHE, 
} ;


/* vacuous declarations and typedef names */

/*     1. Symmetry_t attributes */
typedef enum Symmetry_e     Symmetry_t ;



/* Test the symmetry */
#define Symmetry_IsCylindrical(SYM) \
        (SYM == SYMMETRY_AXIS)
        
#define Symmetry_IsSpherical(SYM) \
        (SYM == SYMMETRY_SPHE)
        
#define Symmetry_IsPlane(SYM) \
        (SYM == SYMMETRY_PLAN)



/* Set the symmetry */
#define Symmetry_SetNoSymmetry(SYM) \
        (SYM = SYMMETRY_NO)
        
#define Symmetry_SetPlaneSymmetry(SYM) \
        (SYM = SYMMETRY_PLAN)
        
#define Symmetry_SetCylindricalSymmetry(SYM) \
        (SYM = SYMMETRY_AXIS)
        
#define Symmetry_SetSphericalSymmetry(SYM) \
        (SYM = SYMMETRY_SPHE)



/* Old notations which I try to eliminate little by little */
#define geom_t    Symmetry_t

#endif
