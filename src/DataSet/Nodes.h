#ifndef NODES_H
#define NODES_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Nodes_s        ; typedef struct Nodes_s        Nodes_t ;


#include "DataFile.h"

extern Nodes_t*  (Nodes_New)                         (const int,const int,const int) ;
extern void      (Nodes_Delete)                      (void*) ;
extern void      (Nodes_CreateMore)                  (Nodes_t*) ;
extern int       (Nodes_ComputeNbOfUnknownFields)    (Nodes_t*) ;
extern void      (Nodes_InitializeObValIndexes)      (Nodes_t*) ;
extern void      (Nodes_SetMatrixRowColumnIndexes)   (Nodes_t*,DataFile_t*) ;
extern void      (Nodes_InitializeMatrixRowColumnIndexes)(Nodes_t*) ;


/* Some constants */
#define Nodes_MaxNbOfMatrices        (10)


#define Nodes_GetNbOfNodes(NODS)          ((NODS)->n_no)
#define Nodes_GetNbOfConnectivities(NODS) ((NODS)->n_con)
#define Nodes_GetNode(NODS)               ((NODS)->no)
#define Nodes_GetNbOfMatrixRows(NODS)     ((NODS)->n_rows)
#define Nodes_GetNbOfMatrixColumns(NODS)  ((NODS)->n_cols)
#define Nodes_GetNbOfDOF(NODS)            ((NODS)->n_dof)
#define Nodes_GetObjectiveValues(NODS)    ((NODS)->obvals)
#define Nodes_GetNbOfMatrices(NODS)       ((NODS)->NbOfMatrices)
#define Nodes_GetPointerToElement(NODS)   ((NODS)->pelement)
#define Nodes_GetBuffers(NODS)            ((NODS)->buffers)



#define Nodes_InitializePointerToElements(NODS) \
        do { \
          int nno = Nodes_GetNbOfNodes(NODS) ; \
          Node_t* node = Nodes_GetNode(NODS) ; \
          Element_t** pel = Node_GetPointerToElement(node) ; \
          int nc = 0 ; \
          int jn ; \
          for(jn = 0 ; jn < nno ; jn++) { \
            Node_t* node_j = node + jn ; \
            Node_GetPointerToElement(node_j) = pel + nc ; \
            nc += Node_GetNbOfElements(node_j) ; \
          } \
        } while(0)




#include "Node.h"
#include "ObVals.h"
#include "Element.h"
#include "Buffers.h"


struct Nodes_s {              /* nodes */
  int NbOfMatrices ;          /* nb of matrices */
  unsigned int* n_rows ;      /* nb of matrix rows */
  unsigned int* n_cols ;      /* nb of matrix columns */
  unsigned int n_no ;         /* nb of nodes */
  unsigned int n_con ;        /* nb of connectivities */
  unsigned int n_dof ;        /* nb of degrees of freedom */
  Node_t* no ;                /* node */
  ObVals_t* obvals ;          /* Objective values */
  Element_tt* pelement ;      /* pointers to elements */
  Buffers_t*  buffers ;       /* Buffers */
} ;

#endif
