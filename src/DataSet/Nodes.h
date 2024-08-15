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


#define Nodes_GetNbOfNodes(NODS)          ((NODS)->NbOfNodes)
#define Nodes_GetNbOfConnectivities(NODS) ((NODS)->NbOfConnectivities)
#define Nodes_GetNode(NODS)               ((NODS)->Node)
#define Nodes_GetNbOfMatrixRows(NODS)     ((NODS)->NbOfMatrixRows)
#define Nodes_GetNbOfMatrixColumns(NODS)  ((NODS)->NbOfMatrixColumns)
#define Nodes_GetNbOfDOF(NODS)            ((NODS)->NbOfDOF)
#define Nodes_GetObjectiveValues(NODS)    ((NODS)->ObjectiveValues)
#define Nodes_GetNbOfMatrices(NODS)       ((NODS)->NbOfMatrices)
#define Nodes_GetPointerToElement(NODS)   ((NODS)->PointerToElement)
#define Nodes_GetBuffers(NODS)            ((NODS)->Buffers)



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


struct Nodes_s {
  Node_t* Node ;
  Element_tt* PointerToElement ;
  ObVals_t* ObjectiveValues ;
  Buffers_t*  Buffers ;
  int NbOfMatrices ;
  unsigned int* NbOfMatrixRows ;
  unsigned int* NbOfMatrixColumns ;
  unsigned int  NbOfNodes ;
  unsigned int  NbOfConnectivities ;
  unsigned int  NbOfDOF ;
} ;

#endif
