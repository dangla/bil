#ifndef FIELDS_H
#define FIELDS_H

/* vacuous declarations and typedef names */

/* class-like structure */
struct Fields_s       ; typedef struct Fields_s       Fields_t ;
/*   1. Fields_t attributes */
struct Field_s          ; typedef struct Field_s          Field_t ;
struct FieldAffine_s    ; typedef struct FieldAffine_s    FieldAffine_t ;
struct FieldGrid_s      ; typedef struct FieldGrid_s      FieldGrid_t ;
struct FieldConstant_s  ; typedef struct FieldConstant_s  FieldConstant_t ;



/* 1. Fields */
#include "Materials.h"
#include "DataFile.h"
#include "Geometry.h"

extern Fields_t* Fields_Create(DataFile_t*,Materials_t*,Geometry_t*) ;

#define Fields_GetNbOfFields(FLDS)    ((FLDS)->n_ch)
#define Fields_GetField(FLDS)         ((FLDS)->ch)




/* 2. Field */

extern double    Field_ComputeValueAtPoint(Field_t*,double*,int) ;

#define Field_MaxLengthOfKeyWord        (30)
#define Field_MaxLengthOfFileName       (200)
#define Field_MaxLengthOfTextLine       (500)

#define Field_GetType(FLD)            ((FLD)->type)
#define Field_GetFieldFormat(FLD)     ((FLD)->store)



/* 3. FieldAffine */

#define FieldAffine_GetValue(FLD)        ((FLD)->v)
#define FieldAffine_GetGradient(FLD)     ((FLD)->g)
#define FieldAffine_GetCoordinate(FLD)   ((FLD)->x)



/* 4. FieldGrid */

#define FieldGrid_GetFileName(FLD)                ((FLD)->name)
#define FieldGrid_GetNbOfPointsAlongX(FLD)        ((FLD)->n_x)
#define FieldGrid_GetNbOfPointsAlongY(FLD)        ((FLD)->n_y)
#define FieldGrid_GetNbOfPointsAlongZ(FLD)        ((FLD)->n_z)
#define FieldGrid_GetCoordinateAlongX(FLD)        ((FLD)->x)
#define FieldGrid_GetCoordinateAlongY(FLD)        ((FLD)->y)
#define FieldGrid_GetCoordinateAlongZ(FLD)        ((FLD)->z)
#define FieldGrid_GetValue(FLD)                   ((FLD)->v)



/* 3. FieldConstant */

#define FieldConstant_GetValue(FLD)      ((FLD)->v)



struct Fields_s {             /* fields */
  unsigned int n_ch ;         /* nb of fields */
  Field_t* ch ;               /* field */
} ;

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
} ;

/* Old notations which I try to eliminate little by little */
#define chmp_t          Field_t
#define chmpaffine_t    FieldAffine_t
#define chmpgrille_t    FieldGrid_t
#define champ(a,b,c)    Field_ComputeValueAtPoint(&(c),(a),(b))

#endif
