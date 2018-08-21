#ifndef CCSFORMAT_H
#define CCSFORMAT_H



/* class-like structure "CCSformat_t" and attributes */

/* vacuous declarations and typedef names */
struct CCSformat_s      ; typedef struct CCSformat_s      CCSformat_t ;



/** The getters */
#define CCSformat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CCSformat_GetNonZeroValue(a)                    ((a)->nzval)
#define CCSformat_GetRowIndexOfTheNonZeroValue(a)       ((a)->rowind)
#define CCSformat_GetFirstNonZeroValueIndexOfColumn(a)  ((a)->colptr)
                                                      




/* complete the structure types by using the typedef */

/* Compressed column storage format (known as Harwell-Boeing sparse matrix format)
 * If a_ij = nzval[k] then rowind[k] = i and colptr[j] <= k < colptr[j + 1] */
struct CCSformat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* rowind ;      /* Row indices of the non zeros */
  unsigned int* colptr ;      /* Index of element in nzval which starts a column */
} ;

#endif
