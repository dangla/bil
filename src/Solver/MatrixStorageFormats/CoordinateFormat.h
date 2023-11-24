#ifndef COORDINATEFORMAT_H
#define COORDINATEFORMAT_H



/* class-like structure "CoordinateFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct CoordinateFormat_s      ; typedef struct CoordinateFormat_s      CoordinateFormat_t ;



#include "Mesh.h"
#include "Options.h"

extern CoordinateFormat_t* (CoordinateFormat_Create)(Mesh_t*,Options_t*,const int) ;
extern void                (CoordinateFormat_Delete)(void*) ;
extern int                 (CoordinateFormat_AssembleElementMatrix)(CoordinateFormat_t*,double*,int*,int*,int,int) ;
extern void                (CoordinateFormat_PrintMatrix)(CoordinateFormat_t*,const int,const char*) ;





/** The getters */
#define CoordinateFormat_GetNbOfNonZeroValues(F)               ((F)->nnz)
#define CoordinateFormat_GetLengthOfArrayValue(F)              ((F)->lvalue)
#define CoordinateFormat_GetLengthOfArrayIndex(F)              ((F)->lindex)
#define CoordinateFormat_GetNonZeroValue(F)                    ((F)->value)
#define CoordinateFormat_GetIndex(F)                           ((F)->index)
//#define CoordinateFormat_GetOptions(F)                         ((F)->options)
                                                      
                                                      
                                                      
                                                      
#define CoordinateFormat_GetColumnIndexOfValue(F) \
        (CoordinateFormat_GetIndex(F) + CoordinateFormat_GetNbOfNonZeroValues(F))


#define CoordinateFormat_GetRowIndexOfValue(F) \
        CoordinateFormat_GetIndex(F)


/* complete the structure types by using the typedef */

/* Coordinate format:
 * a_ij = val[k] ; j = colind[k] ; i = rowind[k]   for  k = 0,..,nnz-1 */
struct CoordinateFormat_s {
  int     nnz ;       /* Nb of non zero values */
  int     lvalue ;    /* Lentgth of array value */
  int     lindex ;    /* Length of array index */
  double* value ;     /* Values */
  int*    index ;     /* Indices */
  int*    colind ;    /* Column indices of the non zeros */
  int*    rowind ;    /* Row indices of the non zeros */
  //Options_t* options ;
} ;

#endif
