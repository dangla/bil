#ifndef RESIDU_H
#define RESIDU_H

/* class-like structures "Residu_t" and attributes */

/* vacuous declarations and typedef names */
struct Residu_s       ; typedef struct Residu_s       Residu_t ;


#include "Element.h"

extern Residu_t*   (Residu_Create)(const int,const int) ;
extern void        (Residu_Delete)(void*) ;
extern void        (Residu_AssembleElementResidu)(Residu_t*,Element_t*,double*) ;



#define Residu_GetResiduIndex(RS)               ((RS)->index)
#define Residu_GetLengthOfRHS(RS)               ((RS)->n)
#define Residu_GetNbOfRHS(RS)                   ((RS)->n_rhs)
#define Residu_GetRHS(RS)                       ((RS)->rhs)
#define Residu_GetSolution(RS)                  ((RS)->sol)


#include "Node.h"

#if 0
#define Residu_AssembleElementResidu(R,EL,RE) \
        do { \
          if(Element_GetMaterial(EL)) { \
            int  nn  = Element_GetNbOfNodes(EL) ; \
            int  neq = Element_GetNbOfEquations(EL) ; \
            int i ; \
            for(i = 0 ; i < nn ; i++) { \
              Node_t* node_i = Element_GetNode(EL,i) ; \
              int j ; \
              for(j = 0 ; j < neq ; j++) { \
                int ij = i*neq + j ; \
                int ii = Element_GetUnknownPosition(EL)[ij] ; \
                if(ii >= 0) { \
                  int k = Node_GetMatrixColumnIndex(node_i)[ii] ; \
                  if(k >= 0) R[k] += RE[ij] ; \
                } \
              } \
            } \
          } \
        } while(0)
#endif
        
        
/* Initialize the residu */
#define Residu_SetValuesToZero(RS) \
        do { \
          unsigned int k ; \
          for(k = 0 ; k < Residu_GetLengthOfRHS(RS) ; k++) { \
            Residu_GetRHS(RS)[k] = 0. ; \
          } \
        } while(0)
        
        

struct Residu_s {             /* Residu */
  unsigned int index ;        /* Matrix index */
  unsigned int    n ;         /* Length of the rhs */
  unsigned int    n_rhs ;     /* Nb of right hand sides */
  double*         rhs ;       /* The values of the rhs */
  double*         sol ;       /* The values of the solutions */
} ;

#endif
