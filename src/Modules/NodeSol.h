#ifndef NODESOL_H
#define NODESOL_H

/* class-like structures "NodesSol_t" and attributes */

/* vacuous declarations and typedef names */
struct NodeSol_s      ; typedef struct NodeSol_s      NodeSol_t ;



extern NodeSol_t* (NodeSol_Create)(const int) ;
extern void       (NodeSol_Delete)(void*) ;
extern NodeSol_t* (NodeSol_GetNodeSolInDistantPast)(NodeSol_t*,unsigned int) ;
extern NodeSol_t* (NodeSol_GetNodeSolInDistantFuture)(NodeSol_t*,unsigned int) ;


#define NodeSol_GetNbOfUnknowns(NS)       ((NS)->nu)
#define NodeSol_GetUnknown(NS)            ((NS)->u)
#define NodeSol_GetPreviousNodeSol(NS)    ((NS)->prev)
#define NodeSol_GetNextNodeSol(NS)        ((NS)->next)






/* Access to unknowns
 * ------------------ */
#define NodeSol_GetCurrentUnknown(NS) \
        NodeSol_GetUnknown(NS)

#define NodeSol_GetPreviousUnknown(NS) \
        NodeSol_GetUnknown(NodeSol_GetPreviousNodeSol(NS))

#define NodeSol_GetNextUnknown(NS) \
        NodeSol_GetUnknown(NodeSol_GetNextNodeSol(NS))
        
#define NodeSol_GetUnknownInDistantPast(NS,i) \
        NodeSol_GetUnknown(NodeSol_GetNodeSolInDistantPast(NS,i))

#define NodeSol_GetUnknownInDistantFuture(NS,i) \
        NodeSol_GetUnknown(NodeSol_GetNodeSolInDistantFuture(NS,i))



struct NodeSol_s {            /* Nodal Solutions */
  unsigned short int nu ;     /* Nb of unknowns */
  double* u ;                 /* Nodal Unknowns */
  NodeSol_t* prev ;           /* Previous Nodal Solutions */
  NodeSol_t* next ;           /* Next Nodal Solutions */
} ;

#endif
