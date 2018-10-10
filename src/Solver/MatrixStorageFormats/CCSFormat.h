#ifndef CCSFORMAT_H
#define CCSFORMAT_H



/* class-like structure "CCSFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct CCSFormat_s      ; typedef struct CCSFormat_s      CCSFormat_t ;



/** The getters */
#define CCSFormat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CCSFormat_GetNonZeroValue(a)                    ((a)->nzval)
#define CCSFormat_GetRowIndexOfNonZeroValue(a)          ((a)->rowind)
#define CCSFormat_GetFirstNonZeroValueIndexOfColumn(a)  ((a)->colptr)
                                                      




/* complete the structure types by using the typedef */

/* Compressed column storage format (known as Harwell-Boeing sparse matrix format)
 * nzval contains the non zero values column-wise
 * rowind contains the row indices
 * colptr points to the columns in both rowind and nzval
 * If a_ij = nzval[k] then rowind[k] = i and colptr[j] <= k < colptr[j + 1] */
struct CCSFormat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* rowind ;      /* Row indices of the non zeros */
  unsigned int* colptr ;      /* Index of element in nzval which starts a column */
} ;

#endif
