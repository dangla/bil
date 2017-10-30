#ifndef NODESSOL_H
#define NODESSOL_H

/* class-like structures "NodesSol_t" and attributes */

/* vacuous declarations and typedef names */
struct NodesSol_s     ; typedef struct NodesSol_s     NodesSol_t ;
struct NodeSol_s      ; typedef struct NodeSol_s      NodeSol_t ;



/* Declaration of Macros, Methods and Structures */

/* 1. NodesSol_t
 * -------------*/
#include "Mesh.h"

extern NodesSol_t*    (NodesSol_Create)(Mesh_t*) ;
extern void           (NodesSol_AllocateMemory)(NodesSol_t*) ;
extern void           (NodesSol_Copy)(NodesSol_t*,NodesSol_t*) ;
 
 
 
#define NodesSol_GetNbOfDOF(NSS)               ((NSS)->NbOfDOF)
#define NodesSol_GetNbOfNodes(NSS)             ((NSS)->NbOfNodes)
#define NodesSol_GetNodeSol(NSS)               ((NSS)->nodesol)





/* Acces to the dof */
#define NodesSol_GetDOF(NSS) \
        NodeSol_GetUnknown(NodesSol_GetNodeSol(NSS))



/* 2. NodeSol_t 
 * ------------*/

extern NodeSol_t* (NodeSol_GetDeepNodeSol)(NodeSol_t*,unsigned int) ;

#define NodeSol_GetNbOfUnknowns(NS)       ((NS)->nu)
#define NodeSol_GetUnknown(NS)            ((NS)->u)
#define NodeSol_GetPreviousNodeSol(NS)    ((NS)->prev)




/* complete the structure types by using the typedef */
struct NodesSol_s {           /* Nodal Solutions */
  unsigned int NbOfNodes ;
  unsigned int NbOfDOF ;      /* Nb of DOF */
  NodeSol_t* nodesol ;
} ;

struct NodeSol_s {            /* Nodal Solutions */
  unsigned short int nu ;     /* Nb of unknowns */
  double* u ;                 /* Nodal Unknowns */
  NodeSol_t* prev ;           /* Previous Nodal Solutions */
} ;

#endif
