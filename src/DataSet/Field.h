#ifndef FIELD_H
#define FIELD_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* vacuous declarations and typedef names */

/* class-like structure */
struct Field_s          ; typedef struct Field_s          Field_t ;
struct FieldAffine_s    ; typedef struct FieldAffine_s    FieldAffine_t ;
struct FieldGrid_s      ; typedef struct FieldGrid_s      FieldGrid_t ;
struct FieldConstant_s  ; typedef struct FieldConstant_s  FieldConstant_t ;


#include "DataFile.h"

extern Field_t*       (Field_New)                  (void) ;
extern void           (Field_Delete)               (void*) ;
extern void           (Field_Scan)                 (Field_t*,DataFile_t*) ;
extern double         (Field_ComputeValueAtPoint)  (Field_t*,double*,int) ;


#define Field_MaxLengthOfKeyWord        (30)
#define Field_MaxLengthOfFileName       (200)
#define Field_MaxLengthOfTextLine       (500)


#define Field_GetType(FLD)            ((FLD)->type)
#define Field_GetFieldFormat(FLD)     ((FLD)->store)



/* 3. FieldAffine */

extern FieldAffine_t* (FieldAffine_Create) (void) ;
extern void           (FieldAffine_Delete) (void*) ;


#define FieldAffine_GetValue(FLD)        ((FLD)->v)
#define FieldAffine_GetGradient(FLD)     ((FLD)->g)
#define FieldAffine_GetCoordinate(FLD)   ((FLD)->x)



/* 4. FieldGrid */

extern FieldGrid_t*   (FieldGrid_Create)(char*) ;
extern void           (FieldGrid_Delete)(void*) ;


#define FieldGrid_GetFileName(FLD)                ((FLD)->name)
#define FieldGrid_GetNbOfPointsAlongX(FLD)        ((FLD)->n_x)
#define FieldGrid_GetNbOfPointsAlongY(FLD)        ((FLD)->n_y)
#define FieldGrid_GetNbOfPointsAlongZ(FLD)        ((FLD)->n_z)
#define FieldGrid_GetCoordinateAlongX(FLD)        ((FLD)->x)
#define FieldGrid_GetCoordinateAlongY(FLD)        ((FLD)->y)
#define FieldGrid_GetCoordinateAlongZ(FLD)        ((FLD)->z)
#define FieldGrid_GetValue(FLD)                   ((FLD)->v)



/* 3. FieldConstant */

#define FieldConstant_GetValue(FLD)               ((FLD)->v)
#define FieldConstant_GetRandomRangeLength(FLD)   ((FLD)->ranlen)



struct Field_s {              /* champ */
  char*   type ;              /* type de champ */
  void*   store ;             /* pointe sur le format du champ */
} ;

struct FieldAffine_s {        /* champ affine */
  double v ;                  /* valeur en A */
  double* g ;                 /* gradient en A */
  double* x ;                 /* coordonnees de A */
} ;

struct FieldGrid_s {          /* champ sur une grille */
  char*  name ;               /* File name which grid is stored in */
  int    n_x ;                /* nb de points sur Ox */
  int    n_y ;                /* nb de points sur Oy */
  int    n_z ;                /* nb de points sur Oz */
  double* x ;                 /* coordonnees sur Ox  */
  double* y ;                 /* coordonnees sur Oy  */
  double* z ;                 /* coordonnees sur Oz  */
  double* v ;                 /* valeurs aux points de la grille */
} ;

struct FieldConstant_s {      /* Constant field */
  double v ;                  /* Value */
  double ranlen ;
} ;

/* Old notations which I try to eliminate little by little */
#define chmp_t          Field_t
#define chmpaffine_t    FieldAffine_t
#define chmpgrille_t    FieldGrid_t
#define champ(a,b,c)    Field_ComputeValueAtPoint(&(c),(a),(b))


#ifdef __CPLUSPLUS
}
#endif
#endif
