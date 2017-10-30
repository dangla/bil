#ifndef SHAPEFCTS_H
#define SHAPEFCTS_H


/* vacuous declarations and typedef names */

/* class-like structure "DataSet_t and attributes */
struct ShapeFcts_s      ; typedef struct ShapeFcts_s      ShapeFcts_t ;

/*   1. ShapeFcts_t attributes */
struct ShapeFct_s       ; typedef struct ShapeFct_s       ShapeFct_t ;


/* Declaration of Macros, Methods and Structures */


/* 1. ShapeFcts_t */
extern ShapeFcts_t*  (ShapeFcts_Create)(void) ;
extern int           (ShapeFcts_FindShapeFct)(ShapeFcts_t*,int,int) ;
extern int           (ShapeFcts_AddShapeFct)(ShapeFcts_t*,int,int) ;


#define ShapeFcts_MaxNbOfShapeFcts             (3)

#define ShapeFcts_GetNbOfShapeFcts(SFS)    ((SFS)->n_sh)
#define ShapeFcts_GetShapeFct(SFS)         ((SFS)->sh)



struct ShapeFcts_s {          /* Shape functions */
  unsigned int n_sh ;         /* Number of shape functions */
  ShapeFct_t*  sh ;           /* Shape function */
} ;



/* 2. ShapeFct_t */
extern void  ShapeFct_ComputeValuesAtPoint(int,int,double*,double*,double*) ;


#define ShapeFct_GetType(SF)             ((SF)->type)
#define ShapeFct_GetNbOfFunctions(SF)    ((SF)->nn)
#define ShapeFct_GetCoordinate(SF)       ((SF)->a)
#define ShapeFct_GetFunction(SF)         ((SF)->h)
#define ShapeFct_GetFunctionGradient(SF) ((SF)->dh)
#define ShapeFct_GetDimension(SF)        ((SF)->dim)



struct ShapeFct_s {           /* Shape function */
  unsigned short int dim ;    /* Sub-dimension (0,1,2,3) */
  unsigned short int nn ;     /* Number of functions */
  double* a ;                 /* Reference coordinates */
  double* h ;                 /* Values of shape functions */
  double* dh ;                /* Values of function gradients */
} ;


#endif
