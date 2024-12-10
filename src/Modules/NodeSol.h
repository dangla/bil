#ifndef NODESOL_H
#define NODESOL_H

#ifdef __CPLUSPLUS
extern "C" {
#endif

/* class-like structures "NodesSol_t" and attributes */

/* vacuous declarations and typedef names */
struct NodeSol_s      ; typedef struct NodeSol_s      NodeSol_t ;



extern NodeSol_t* (NodeSol_Create)(const int) ;
extern void       (NodeSol_Delete)(void*) ;
extern void       (NodeSol_Copy)(NodeSol_t*,NodeSol_t*) ;


#define NodeSol_GetNbOfUnknowns(NS)       ((NS)->nu)
#define NodeSol_GetUnknown(NS)            ((NS)->u)
//#define NodeSol_GetPreviousNodeSol(NS)    ((NS)->prev)
//#define NodeSol_GetNextNodeSol(NS)        ((NS)->next)




struct NodeSol_s {            /* Nodal Solutions */
  unsigned int nu ;     /* Nb of unknowns */
  double* u ;                 /* Nodal Unknowns */
  //NodeSol_t* prev ;           /* Previous Nodal Solutions */
  //NodeSol_t* next ;           /* Next Nodal Solutions */
} ;


#ifdef __CPLUSPLUS
}
#endif
#endif
