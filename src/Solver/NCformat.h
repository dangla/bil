#ifndef NCFORMAT_H
#define NCFORMAT_H


/* class-like structure "NCformat_t" and attributes */

/* vacuous declarations and typedef names */
//struct NCformat_s      ; typedef struct NCformat_s      NCformat_t ;
#define NCformat_t        NCformat     /* Mimic the NCformat struct defined in SuperLu */


#include "Mesh.h"

#ifdef SLU_DIR
  #define val(x) #x
  #define xval(x) val(x)
  #include xval(SLU_DIR/SRC/dsp_defs.h)
  #undef val
  #undef xval

  extern NCformat_t* (NCformat_Create)(Mesh_t*) ;
  extern void        (NCformat_Delete)(void*) ;
  extern void (NCformat_AssembleElementMatrix)(NCformat_t*,double*,int*,int*,int,int*,int) ;
  extern void (NCformat_PrintMatrix)(NCformat_t*,unsigned int,const char*) ;
#endif



/** The getters */
#define NCformat_GetNbOfNonZeroValues(a)               ((a)->nnz)
#define NCformat_GetNonZeroValue(a)                    ((a)->nzval)
#define NCformat_GetRowIndexOfTheNonZeroValue(a)       ((a)->rowind)
#define NCformat_GetFirstNonZeroValueIndexOfColumn(a)  ((a)->colptr)

/*
#define NCformat_NbOfNonZeroValuesInColumn(a,j)        ((a)->colptr[j + 1] - (a)->colptr[j])
#define NCformat_RowIndexStartingColumn(a,j)           ((a)->rowind[(a)->colptr[j]])
#define NCformat_UpperColumn(a,j)                      ((a)->nzval + (a)->colptr[j + 1] - (j))
*/


/* complete the structure types by using the typedef */

/* Same as CCSformat
 * Compressed column storage format (known as Harwell-Boeing sparse matrix format)
 * If a_ij = nzval[k] then rowind[k] = i and colptr[j] <= k < colptr[j + 1] */
struct NCformat_s {
  unsigned int    nnz ;       /* nb of non zero values */
  double* nzval ;             /* Non zero values */
  unsigned int* rowind ;      /* Row indices of the non zeros */
  unsigned int* colptr ;      /* Index of element in nzval which starts a column */
} ;


#endif
