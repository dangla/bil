#ifndef MATRIXSTORAGEFORMAT_H
#define MATRIXSTORAGEFORMAT_H


enum MatrixFormat_e {         /* format of the matrix to be stored */
  LDUSKL,                     /* LDU Skyline format */
  HBSLU,                      /* Harwell-Boeing SuperLU format */
  CCS = HBSLU,                /* Compressed column storage format */
  CRS                         /* Compressed row storage format */
} ;


/* class-like structure "MatrixStorageFormat_t" and attributes */

/* vacuous declarations and typedef names */
typedef enum MatrixFormat_e MatrixFormat_t ;
struct LDUSKLformat_s   ; typedef struct LDUSKLformat_s   LDUSKLformat_t ;
struct CRSformat_s      ; typedef struct CRSformat_s      CRSformat_t ;
struct CCSformat_s      ; typedef struct CCSformat_s      CCSformat_t ;
/* struct SuperMatrix_s    ; typedef struct SuperMatrix_s    SuperMatrix_t ; */
#define SuperMatrix_t     SuperMatrix  /* Mimic the SuperMatrix struct defined in SuperLu */
#define NCformat_t        NCformat     /* Mimic the NCformat struct defined in SuperLu */


#include "Mesh.h"

#ifdef SLU_DIR
  #define val(x) #x
  #define xval(x) val(x)
  #include xval(SLU_DIR/SRC/dsp_defs.h)
  #undef val
  #undef xval
#endif


extern LDUSKLformat_t* (LDUSKLformat_Create)(Mesh_t*) ;
extern void            (LDUSKLformat_Delete)(void*) ;
extern void LDUSKLformat_AssembleElementMatrix(LDUSKLformat_t*,double*,int*,int*,int) ;
extern void LDUSKLformat_PrintMatrix(LDUSKLformat_t*,unsigned int,const char*) ;

#ifdef SLU_DIR
extern SuperMatrix_t* (SuperMatrix_Create)(Mesh_t*) ;
extern void           (SuperMatrix_Delete)(void*) ;

extern NCformat_t* (NCformat_Create)(Mesh_t*) ;
extern void        (NCformat_Delete)(void*) ;
extern void NCformat_AssembleElementMatrix(NCformat_t*,double*,int*,int*,int,int*,int) ;
extern void NCformat_PrintMatrix(NCformat_t*,unsigned int,const char*) ;
#endif



/** Macros for LDUSKLformat ------------------------------------------*/
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
                                                      



/** Macros for CCSformat --------------------------------------------*/
#define CCSformat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CCSformat_GetNonZeroValue(a)                    ((a)->nzval)
#define CCSformat_GetRowIndexOfTheNonZeroValue(a)       ((a)->rowind)
#define CCSformat_GetFirstNonZeroValueIndexOfColumn(a)  ((a)->colptr)
                                                      



/** Macros for NCformat (same as CCSformat)--------------------------*/
#define NCformat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define NCformat_GetNonZeroValue(a)                    ((a)->nzval)
#define NCformat_GetRowIndexOfTheNonZeroValue(a)       ((a)->rowind)
#define NCformat_GetFirstNonZeroValueIndexOfColumn(a)  ((a)->colptr)

/*
#define NCformat_NbOfNonZeroValuesInColumn(a,j)        ((a)->colptr[j + 1] - (a)->colptr[j])
#define NCformat_RowIndexStartingColumn(a,j)           ((a)->rowind[(a)->colptr[j]])
#define NCformat_UpperColumn(a,j)                      ((a)->nzval + (a)->colptr[j + 1] - (j))
*/



/** Macros for SuperMatrix -------------------------------------------*/
#define SuperMatrix_GetNbOfRows(a)                     ((a)->nrow)
#define SuperMatrix_GetNbOfColumns(a)                  ((a)->ncol)
#define SuperMatrix_GetStorage(a)                      ((a)->Store)
#define SuperMatrix_GetStorageType(a)                  ((a)->Stype)
#define SuperMatrix_GetDataType(a)                     ((a)->Dtype)
#define SuperMatrix_GetMatrixType(a)                   ((a)->Mtype)



/* complete the structure types by using the typedef */

/* LDU Skyline format */
struct LDUSKLformat_s {       /* LDU Skyline storage format */
  unsigned int    nnz ;       /* nb of non zero values */
  double* d ;                 /* diagonal matrix values */
  double** l ;                /* pointer to strictly lower triangular matrix values */
  double** u ;                /* pointer to strictly upper triangular matrix values */
} ;



/* Skyline format */
struct SKLformat_s {          /* Skyline storage format */
  unsigned int    nnz ;       /* nb of non zero values */
  double* d ;                 /* diagonal matrix values */
  double* l ;                 /* strictly lower triangular matrix values */
  double* u ;                 /* strictly upper triangular matrix values */
  unsigned int* colptr ;      /* Index of element in u which starts a column */
  unsigned int* rowptr ;      /* Index of element in l which starts a row */
} ;



/* Compressed row storage format
 * If a_ij = nzval[k] then colind[k] = j and rowptr[i] <= k < rowptr[i + 1] */
struct CRSformat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* colind ;      /* Column indices of the non zeros */
  unsigned int* rowptr ;      /* Index of element in nzval which starts a row */
} ;



/* Compressed column storage format (known as Harwell-Boeing sparse matrix format)
 * If a_ij = nzval[k] then rowind[k] = i and colptr[j] <= k < colptr[j + 1] */
struct CCSformat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* rowind ;      /* Row indices of the non zeros */
  unsigned int* colptr ;      /* Index of element in nzval which starts a column */
} ;



/* Copy of the SuperMatrix struc as defined in SuperLU */
#ifdef SLU_DIR
struct SuperMatrix_s {
  Stype_t Stype;          /* Storage format (SLU_[NC,NR,SC,SR,NCP,DN]) */
  Dtype_t Dtype;          /* Data type (SLU_[S,D,C,Z])*/
  Mtype_t Mtype;          /* Mathematical property of the matrix (SLU_[GE,..]) */
  int    nrow;            /* number of rows */
  int    ncol;            /* number of columns */
  void*  Store;           /* pointer to the actual storage of the matrix */
} ;
#endif


#endif
