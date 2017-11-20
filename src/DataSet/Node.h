#ifndef NODE_H
#define NODE_H


/* vacuous declarations and typedef names */

/* class-like structure */
struct Node_s         ; typedef struct Node_s         Node_t ;



/* Accessors */
#define Node_GetNodeIndex(NOD)            ((NOD)->index)
#define Node_GetCoordinate(NOD)           ((NOD)->x)
#define Node_GetNbOfEquations(NOD)        ((NOD)->neq)
#define Node_GetNbOfUnknowns(NOD)         ((NOD)->nin)
#define Node_GetNameOfEquation(NOD)       ((NOD)->eqn)
#define Node_GetNameOfUnknown(NOD)        ((NOD)->inc)
#define Node_GetObValIndex(NOD)           ((NOD)->i_obj)
#define Node_GetMatrixColumnIndex(NOD)    ((NOD)->colindex)
#define Node_GetMatrixRowIndex(NOD)       ((NOD)->rowindex)
#define Node_GetNodeSol(NOD)              ((NOD)->sol)



/* Access to NodeSol */
#define Node_GetPreviousNodeSol(NOD) \
        NodeSol_GetPreviousNodeSol(Node_GetNodeSol(NOD))

#define Node_GetDeepNodeSol(NOD,i) \
        NodeSol_GetDeepNodeSol(Node_GetNodeSol(NOD),i)



/* Access to unknowns*/
#define Node_GetCurrentUnknown(NOD) \
        NodeSol_GetUnknown(Node_GetNodeSol(NOD))

#define Node_GetPreviousUnknown(NOD) \
        NodeSol_GetUnknown(Node_GetPreviousNodeSol(NOD))

#define Node_GetDeepUnknown(NOD,i) \
        NodeSol_GetUnknown(Node_GetDeepNodeSol(NOD,i))




#include "NodeSol.h"


struct Node_s {               /* noeud */
  unsigned int index ;        /* node index */
  double* x ;                 /* coordonnees */
  unsigned short int neq ;    /* nombre d'equations au noeud */
  unsigned short int nin ;    /* nombre d'inconnues au noeud */
  char**   inc ;              /* nom des inconnues */
  char**   eqn ;              /* nom des equations */
  unsigned short int* i_obj ; /* indices des inconnues dans obj */
  int*    colindex ;          /* column index (unknowns) */
  int*    rowindex ;          /* row index (equations) */
  NodeSol_t* sol ;            /* Nodal solution */
} ;

#endif
