#ifndef RESIDU_H
#define RESIDU_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "Residu_t" and attributes */

/* vacuous declarations and typedef names */
struct Residu_s       ; typedef struct Residu_s       Residu_t ;


#include "Mesh.h"
#include "Options.h"
#include "Element.h"

extern Residu_t*   (Residu_Create)(Mesh_t*,Options_t*,const int,const int) ;
extern void        (Residu_Delete)(void*) ;
extern int         (Residu_AssembleElementResidu)(Residu_t*,Element_t*,double*) ;
extern void        (Residu_PrintResidu)(Residu_t*,const char*) ;
extern void        (Residu_SetValuesToZero)(Residu_t*) ;



#define Residu_GetResiduIndex(RS)               ((RS)->index)
#define Residu_GetVectorStorageFormat(RS)       ((RS)->vectorstorageformat)
#define Residu_GetLengthOfRHS(RS)               ((RS)->n)
#define Residu_GetNbOfRHS(RS)                   ((RS)->n_rhs)
#define Residu_GetNbOfSolutions(RS)             ((RS)->n_rhs)
#define Residu_GetRHS(RS)                       ((RS)->rhs)
#define Residu_GetSolution(RS)                  ((RS)->sol)
#define Residu_GetStoragOfRHS(RS)               ((RS)->storageofrhs)
#define Residu_GetStoragOfSolution(RS)          ((RS)->storageofsol)


#define Residu_GetOptions(RS) \
        VectorStorageFormat_GetOptions(Residu_GetVectorStorageFormat(RS))

#include "Node.h"
        
        
/* Initialize the residu */
#if 0
#define Residu_SetValuesToZero(RS) \
        do { \
          unsigned int Residu_n = Residu_GetNbOfRHS(RS)*Residu_GetLengthOfRHS(RS) ; \
          double* Residu_d = (double*) Residu_GetRHS(RS) ; \
          unsigned int Residu_k ; \
          for(Residu_k = 0 ; Residu_k < Residu_n ; Residu_k++) { \
            Residu_d[Residu_k] = 0. ; \
          } \
        } while(0)
#endif



#include "VectorStorageFormat.h"

#define Residu_GetStorageFormatType(RS) \
        VectorStorageFormat_GetType(Residu_GetVectorStorageFormat(RS))

#define Residu_StorageFormatIs(RS,KEY) \
        VectorStorageFormat_Is(Residu_GetVectorStorageFormat(RS),KEY)


struct Residu_s {             /* Residu */
  unsigned int index ;        /* Residu index */
  VectorStorageFormat_t* vectorstorageformat ; /* Storage format of vectors */
  unsigned int    n ;         /* Length of the rhs */
  unsigned int    n_rhs ;     /* Nb of right hand sides */
  double*         rhs ;       /* The values of the rhs */
  double*         sol ;       /* The values of the solutions */
  void*           storageofrhs ;
  void*           storageofsol ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
