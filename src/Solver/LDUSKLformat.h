#ifndef LDUSKLFORMAT_H
#define LDUSKLFORMAT_H



/* class-like structure "LDUSKLformat_t" and attributes */

/* vacuous declarations and typedef names */
struct LDUSKLformat_s   ; typedef struct LDUSKLformat_s   LDUSKLformat_t ;


#include "Mesh.h"


extern LDUSKLformat_t* (LDUSKLformat_Create)(Mesh_t*) ;
extern void            (LDUSKLformat_Delete)(void*) ;
extern void LDUSKLformat_AssembleElementMatrix(LDUSKLformat_t*,double*,int*,int*,int) ;
extern void LDUSKLformat_PrintMatrix(LDUSKLformat_t*,unsigned int,const char*) ;




/** The getters */
#define LDUSKLformat_GetNbOfNonZeroValues(a)         ((a)->nnz)
#define LDUSKLformat_GetNonZeroValue(a)              ((a)->d)
#define LDUSKLformat_GetDiagonal(a)                  ((a)->d)
#define LDUSKLformat_GetPointerToLowerRow(a)         ((a)->l)
#define LDUSKLformat_GetPointerToUpperColumn(a)      ((a)->u)



#define LDUSKLformat_GetUpperColumn(a,j) \
        (LDUSKLformat_GetPointerToUpperColumn(a)[j] - (j))

#define LDUSKLformat_GetLowerRow(a,i) \
        (LDUSKLformat_GetPointerToLowerRow(a)[i] - (i))

#define LDUSKLformat_HeightOfUpperColumn(a,j) \
        ((j > 0) ? (LDUSKLformat_GetPointerToUpperColumn(a)[j] - \
                    LDUSKLformat_GetPointerToUpperColumn(a)[j - 1]) : 0)

#define LDUSKLformat_LengthOfLowerRow(a,i) \
        ((i > 0) ? (LDUSKLformat_GetPointerToLowerRow(a)[i] - \
                    LDUSKLformat_GetPointerToLowerRow(a)[i - 1]) : 0)

/*
#define LDUSKLformat_UpperValue(a,i,j)     (LDUSKLformat_GetUpperColumn(a,j)[i])
#define LDUSKLformat_LowerValue(a,i,j)     (LDUSKLformat_GetLowerRow(a,i)[j])
#define LDUSKLformat_DiagonalValue(a,i)    (LDUSKLformat_GetDiagonal(a)[i])
#define LDUSKLformat_Value(a,i,j)          ((i < j) ? LDUSKLformat_UpperValue(a,i,j) : \
                                           ((i > j) ? LDUSKLformat_LowerValue(a,i,j) : \
                                                      LDUSKLformat_DiagonalValue(a)[i]))
*/

#define LDUSKLformat_ColumnIndexStartingRow(a,i) \
        ((i > 0) ? (i) - LDUSKLformat_LengthOfLowerRow(a,i) : 0)

#define LDUSKLformat_RowIndexStartingColumn(a,j) \
        ((j > 0) ? (j) - LDUSKLformat_HeightOfUpperColumn(a,j) : 0)
                                                      



/* complete the structure types by using the typedef */

/* LDU Skyline format */
struct LDUSKLformat_s {       /* LDU Skyline storage format */
  unsigned int    nnz ;       /* nb of non zero values */
  double* d ;                 /* diagonal matrix values */
  double** l ;                /* pointer to strictly lower triangular matrix values */
  double** u ;                /* pointer to strictly upper triangular matrix values */
} ;


#endif
