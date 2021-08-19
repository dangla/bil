#ifndef NODESOL_H
#define NODESOL_H

/* class-like structures "NodesSol_t" and attributes */

/* vacuous declarations and typedef names */
struct NodeSol_s      ; typedef struct NodeSol_s      NodeSol_t ;



extern NodeSol_t* (NodeSol_Create)(const int) ;
extern NodeSol_t* (NodeSol_GetNodeSolInDistantPast)(NodeSol_t*,unsigned int) ;
extern NodeSol_t* (NodeSol_GetNodeSolInDistantFuture)(NodeSol_t*,unsigned int) ;


#define NodeSol_GetNbOfUnknowns(NS)       ((NS)->nu)
#define NodeSol_GetUnknown(NS)            ((NS)->u)
#define NodeSol_GetPreviousNodeSol(NS)    ((NS)->prev)
#define NodeSol_GetNextNodeSol(NS)        ((NS)->next)



struct NodeSol_s {            /* Nodal Solutions */
  unsigned short int nu ;     /* Nb of unknowns */
  double* u ;                 /* Nodal Unknowns */
  NodeSol_t* prev ;           /* Previous Nodal Solutions */
  NodeSol_t* next ;           /* Next Nodal Solutions */
} ;

#endif
