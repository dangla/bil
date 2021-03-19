#ifndef NODES_H
#define NODES_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Nodes_s        ; typedef struct Nodes_s        Nodes_t ;


#include "DataFile.h"

extern Nodes_t*  (Nodes_New)                         (const int,const int,const int) ;
extern void      (Nodes_Delete)                      (void*) ;
extern void      (Nodes_CreateMore)                  (Nodes_t*) ;
extern void      (Nodes_DeleteMore)                  (void*) ;
extern int       (Nodes_ComputeNbOfUnknownFields)    (Nodes_t*) ;
extern void      (Nodes_InitializeObValIndexes)      (Nodes_t*) ;
extern void      (Nodes_SetMatrixRowColumnIndexes)   (Nodes_t*,DataFile_t*) ;
extern void      (Nodes_InitializeMatrixRowColumnIndexes)(Nodes_t*) ;



#define Nodes_GetNbOfNodes(NODS)          ((NODS)->n_no)
#define Nodes_GetNbOfConnectivities(NODS) ((NODS)->n_con)
#define Nodes_GetNode(NODS)               ((NODS)->no)
#define Nodes_GetNbOfMatrixRows(NODS)     ((NODS)->n_rows)
#define Nodes_GetNbOfDOF(NODS)            ((NODS)->n_dof)
#define Nodes_GetObjectiveValues(NODS)    ((NODS)->obvals)


/* Synonyms */
#define Nodes_GetNbOfMatrixColumns \
        Nodes_GetNbOfMatrixRows




#include "Node.h"
#include "ObVals.h"


struct Nodes_s {              /* nodes */
  unsigned int n_rows ;       /* nb of matrix rows */
  unsigned int n_no ;         /* nb of nodes */
  unsigned int n_con ;        /* nb of connectivities */
  unsigned int n_dof ;        /* nb of degrees of freedom */
  Node_t* no ;                /* node */
  ObVals_t* obvals ;          /* Objective values */
} ;

#endif
