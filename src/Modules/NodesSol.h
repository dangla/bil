#ifndef NODESSOL_H
#define NODESSOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "NodesSol_t" and attributes */

/* vacuous declarations and typedef names */
struct NodesSol_s     ; typedef struct NodesSol_s     NodesSol_t ;



#include "Mesh.h"

extern NodesSol_t*    (NodesSol_Create)(Mesh_t*) ;
extern void           (NodesSol_Delete)(void*) ;
extern void           (NodesSol_Copy)(NodesSol_t*,NodesSol_t*) ;
extern void           (NodesSol_CopySelectedSequentialUnknowns)(NodesSol_t*,NodesSol_t*,const int) ;
 
 
 
#define NodesSol_GetNodes(NSS)                 ((NSS)->nodes)
//#define NodesSol_GetNbOfDOF(NSS)               ((NSS)->NbOfDOF)
#define NodesSol_GetNbOfNodes(NSS)             ((NSS)->NbOfNodes)
#define NodesSol_GetNodeSol(NSS)               ((NSS)->nodesol)
//#define NodesSol_GetNodalValue(NSS)            ((NSS)->NodalValue)




#include "Nodes.h"
#include "NodeSol.h"


/* complete the structure types by using the typedef */
struct NodesSol_s {           /* Nodal Solutions */
  Nodes_t* nodes ;
  unsigned int NbOfNodes ;
  //unsigned int NbOfDOF ;      /* Nb of DOF */
  //double*      NodalValue ;
  NodeSol_t* nodesol ;
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
