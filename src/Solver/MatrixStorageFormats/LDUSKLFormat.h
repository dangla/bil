#ifndef LDUSKLFORMAT_H
#define LDUSKLFORMAT_H



/* class-like structure "LDUSKLFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct LDUSKLFormat_s   ; typedef struct LDUSKLFormat_s   LDUSKLFormat_t ;


#include "Mesh.h"


//extern LDUSKLFormat_t* (LDUSKLFormat_Create)(Mesh_t*) ;
extern LDUSKLFormat_t* (LDUSKLFormat_Create)(Mesh_t*,const int) ;
extern void            (LDUSKLFormat_Delete)(void*) ;
extern void LDUSKLFormat_AssembleElementMatrix(LDUSKLFormat_t*,double*,int*,int*,int) ;
extern void LDUSKLFormat_PrintMatrix(LDUSKLFormat_t*,unsigned int,const char*) ;




/** The getters */
#define LDUSKLFormat_GetNbOfNonZeroValues(a)         ((a)->nnz)
#define LDUSKLFormat_GetNonZeroValue(a)              ((a)->d)
#define LDUSKLFormat_GetDiagonal(a)                  ((a)->d)
#define LDUSKLFormat_GetPointerToLowerRow(a)         ((a)->l)
#define LDUSKLFormat_GetPointerToUpperColumn(a)      ((a)->u)



#define LDUSKLFormat_GetUpperColumn(a,j) \
        (LDUSKLFormat_GetPointerToUpperColumn(a)[j] - (j))

#define LDUSKLFormat_GetLowerRow(a,i) \
        (LDUSKLFormat_GetPointerToLowerRow(a)[i] - (i))

#define LDUSKLFormat_HeightOfUpperColumn(a,j) \
        ((j > 0) ? (LDUSKLFormat_GetPointerToUpperColumn(a)[j] - \
                    LDUSKLFormat_GetPointerToUpperColumn(a)[j - 1]) : 0)

#define LDUSKLFormat_LengthOfLowerRow(a,i) \
        ((i > 0) ? (LDUSKLFormat_GetPointerToLowerRow(a)[i] - \
                    LDUSKLFormat_GetPointerToLowerRow(a)[i - 1]) : 0)


/* Components of the matrix */
#define LDUSKLFormat_UpperValue(a,i,j) \
        ((i >= LDUSKLFormat_RowIndexStartingColumn(a,j) ? LDUSKLFormat_GetUpperColumn(a,j)[i] : 0)
       
#define LDUSKLFormat_LowerValue(a,i,j) \
        ((j >= LDUSKLFormat_ColumnIndexStartingRow(a,i)) ? LDUSKLFormat_GetLowerRow(a,i)[j] : 0)
        
#define LDUSKLFormat_DiagonalValue(a,i) \
        (LDUSKLFormat_GetDiagonal(a)[i])
         
#define LDUSKLFormat_Value(a,i,j) \
        ((i < j) ? LDUSKLFormat_UpperValue(a,i,j) : \
        ((i > j) ? LDUSKLFormat_LowerValue(a,i,j) : \
        LDUSKLFormat_DiagonalValue(a,i)))


/* Starting indexes */
#define LDUSKLFormat_ColumnIndexStartingRow(a,i) \
        ((i > 0) ? (i) - LDUSKLFormat_LengthOfLowerRow(a,i) : 0)

#define LDUSKLFormat_RowIndexStartingColumn(a,j) \
        ((j > 0) ? (j) - LDUSKLFormat_HeightOfUpperColumn(a,j) : 0)
                                                      



/* complete the structure types by using the typedef */

/* LDU Skyline format */
struct LDUSKLFormat_s {       /* LDU Skyline storage format */
  unsigned int    nnz ;       /* Nb of non zero values */
  double* d ;                 /* Diagonal matrix values */
  double** l ;                /* Pointer to strictly lower triangular matrix values */
  double** u ;                /* Pointer to strictly upper triangular matrix values */
} ;


#endif
