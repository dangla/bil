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
#define Residu_GetNbOfSolutions(RS)             ((RS)->n_rhs)
#define Residu_GetRHS(RS)                       ((RS)->rhs)
#define Residu_GetSolution(RS)                  ((RS)->sol)


#include "Node.h"

#if 0
#define Residu_AssembleElementResidu(R,EL,RE) \
        do { \
          if(Element_GetMaterial(EL)) { \
            int Residu_nn  = Element_GetNbOfNodes(EL) ; \
            int Residu_neq = Element_GetNbOfEquations(EL) ; \
            int Residu_i ; \
            for(Residu_i = 0 ; Residu_i < Residu_nn ; Residu_i++) { \
              Node_t* node_i = Element_GetNode(EL,Residu_i) ; \
              int Residu_j ; \
              for(Residu_j = 0 ; Residu_j < Residu_neq ; Residu_j++) { \
                int Residu_ij = Residu_i*Residu_neq + Residu_j ; \
                int Residu_ii = Element_GetUnknownPosition(EL)[Residu_ij] ; \
                if(Residu_ii >= 0) { \
                  int Residu_k = Node_GetMatrixColumnIndex(node_i)[Residu_ii] ; \
                  if(Residu_k >= 0) R[Residu_k] += RE[Residu_ij] ; \
                } \
              } \
            } \
          } \
        } while(0)
#endif
        
        
/* Initialize the residu */
#define Residu_SetValuesToZero(RS) \
        do { \
          unsigned int Residu_n = Residu_GetNbOfRHS(RS)*Residu_GetLengthOfRHS(RS) ; \
          unsigned int Residu_k ; \
          for(Residu_k = 0 ; Residu_k < Residu_n ; Residu_k++) { \
            Residu_GetRHS(RS)[Residu_k] = 0. ; \
          } \
        } while(0)
        
        

struct Residu_s {             /* Residu */
  unsigned int index ;        /* Residu index */
  unsigned int    n ;         /* Length of the rhs */
  unsigned int    n_rhs ;     /* Nb of right hand sides */
  double*         rhs ;       /* The values of the rhs */
  double*         sol ;       /* The values of the solutions */
} ;

#endif
