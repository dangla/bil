#ifndef CRSFORMAT_H
#define CRSFORMAT_H



/* class-like structure "CRSFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct CRSFormat_s      ; typedef struct CRSFormat_s      CRSFormat_t ;

#include "Mesh.h"



/** The getters */
#define CRSFormat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define CRSFormat_GetNonZeroValue(a)                    ((a)->nzval)
#define CRSFormat_GetColumnIndexOfNonZeroValue(a)       ((a)->colind)
#define CRSFormat_GetFirstNonZeroValueIndexOfRow(a)     ((a)->rowptr)
                                                      


/* complete the structure types by using the typedef */

/* Compressed row storage format
 * If a_ij = nzval[k] then colind[k] = j and rowptr[i] <= k < rowptr[i + 1] */
struct CRSFormat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* colind ;      /* Column indices of the non zeros */
  unsigned int* rowptr ;      /* Index of element in nzval which starts a row */
} ;

#endif
