#ifndef SKLFORMAT_H
#define SKLFORMAT_H


/* class-like structure "SKLformat_t" and attributes */

/* vacuous declarations and typedef names */
struct SKLformat_s   ; typedef struct SKLformat_s   SKLformat_t ;



/** The getters */
#define SKLformat_GetNbOfNonZeroValues(a)                    ((a)->nnz)
#define SKLformat_GetNonZeroValue(a)                         ((a)->d)
#define SKLformat_GetDiagonal(a)                             ((a)->d)
#define SKLformat_GetPointerToLowerRow(a)                    ((a)->l)
#define SKLformat_GetPointerToUpperColumn(a)                 ((a)->u)
#define SKLformat_GetFirstNonZeroValueIndexOfLowerRow(a)     ((a)->rowptr)
#define SKLformat_GetFirstNonZeroValueIndexOfUpperColumn(a)  ((a)->colptr)



/* complete the structure types by using the typedef */

/* Skyline format */
struct SKLformat_s {          /* Skyline storage format */
  unsigned int    nnz ;       /* nb of non zero values */
  double* d ;                 /* diagonal matrix values */
  double* l ;                 /* strictly lower triangular matrix values */
  double* u ;                 /* strictly upper triangular matrix values */
  unsigned int* colptr ;      /* Index of element in u which starts a column */
  unsigned int* rowptr ;      /* Index of element in l which starts a row */
} ;

#endif
