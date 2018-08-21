#ifndef CRSFORMAT_H
#define CRSFORMAT_H



/* class-like structure "CRSformat_t" and attributes */

/* vacuous declarations and typedef names */
struct CRSformat_s      ; typedef struct CRSformat_s      CRSformat_t ;

#include "Mesh.h"



/** The getters */
#define CRSformat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CRSformat_GetNonZeroValue(a)                    ((a)->nzval)
#define CRSformat_GetColumnIndexOfTheNonZeroValue(a)    ((a)->colind)
#define CRSformat_GetFirstNonZeroValueIndexOfRow(a)     ((a)->rowptr)
                                                      


/* complete the structure types by using the typedef */

/* Compressed row storage format
 * If a_ij = nzval[k] then colind[k] = j and rowptr[i] <= k < rowptr[i + 1] */
struct CRSformat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* colind ;      /* Column indices of the non zeros */
  unsigned int* rowptr ;      /* Index of element in nzval which starts a row */
} ;

#endif
