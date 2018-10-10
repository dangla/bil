#ifndef SKLFORMAT_H
#define SKLFORMAT_H


/* class-like structure "SKLFormat_t" and attributes */

/* vacuous declarations and typedef names */
struct SKLFormat_s   ; typedef struct SKLFormat_s   SKLFormat_t ;



/** The getters */
#define SKLFormat_GetNbOfNonZeroValues(a)                    ((a)->nnz)
#define SKLFormat_GetNonZeroValue(a)                         ((a)->nzval)
#define SKLFormat_GetPointerToLowerRow(a)                    ((a)->l)
#define SKLFormat_GetPointerToUpperColumn(a)                 ((a)->u)
#define SKLFormat_GetFirstNonZeroValueIndexOfLowerRow(a)     ((a)->rowptr)
#define SKLFormat_GetFirstNonZeroValueIndexOfUpperColumn(a)  ((a)->colptr)



/* complete the structure types by using the typedef */

/* Skyline format */
struct SKLFormat_s {          /* Skyline storage format */
  unsigned int    nnz ;       /* Nb of non zero values */
  double* nzval ;             /* Non zero values */
  double* l ;                 /* Strictly lower triangular matrix values */
  double* u ;                 /* Upper triangular matrix values including diagonal */
  unsigned int* colptr ;      /* Index of element in u which starts a column */
  unsigned int* rowptr ;      /* Index of element in l which starts a row */
} ;

#endif
